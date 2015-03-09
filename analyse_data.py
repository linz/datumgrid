#!/usr/bin/python

import sys
import os
import os.path
import pandas as pd
import numpy as np
import argparse
from subprocess import call

parser=argparse.ArgumentParser(description='Iterative analysis of deformation field')
parser.add_argument('input_file',help='Input deformation file')
parser.add_argument('-c','--calcdir',help='Directory for output files')
parser.add_argument('-f','--force-initial-reweighting',action='store_true',help='Force calculation of initial rejections')
parser.add_argument('-i','--ignore-initial-reweighting',action='store_true',help='Ignore initial reweighting in iterative reweight solution')
args=parser.parse_args()

# Input file is a CSV file with columns
# code,lon,lat,de,den,order

def_input_file=args.input_file

calcdir=args.calcdir or os.path.join(os.path.dirname(def_input_file),'calcs')

if not os.path.exists( calcdir ):
    os.makedirs(calcdir)

calcdir = calcdir+'/'

# datafile='model_data.csv'
summary_file=calcdir+'reweight_summary.csv'
analysis_log=calcdir+'analysis.log'
if os.path.exists(analysis_log):
    os.remove(analysis_log)

# Reference value for distortion error (datumgrid parameter)

base_distortion_error=25.0

# Phase 1 parameters
#
# Initial downweigthing of misfitting stations

max_iterations=10  
grid_spacing_initial=0.01 # Initial grid latitude spacing (longitude=latitude*1.5)
grid_spacing_final=0.002  # Final grid latitude spacing (longitude=latitude*1.5)
percent_used=5.0          # Percentage of standardised residuals considered
min_std_res=5.0           # If all values are less than this then don't do anything.
influence_radius=2000.0   # Points within this range of point already downgraded 
                          # in current iteration are not considered.
downgrade_factor=0.25     # Amount by which rejected station is downgraded

# Phase 2 parameters
#
# Iterative reweighting of data .. details below

reweight_iterations=10    # Number of reweighting iterations
min_coord_change=0.01     # Tolerable coordinate change to end iterations
min_residual=0.05         # Minumum residual - use normal weighting below this, 
                          # otherwise reweight to 1-norm to avoid too much influence
max_residual=0.2          # If greater than this then downweight further 


def run_datumgrid( data_file, result_file, point_error=0.05, grid_spacing=0.01, proximity=5000, distortion_error=base_distortion_error,
                downgrades=None, calcstdres=False ):
    dgcmd='''
    # datumgrid control file

    data_file {data_file}
    grid_spacing {grid_spacing}
    coordinate_to_metres 73000 100000
    required_point_proximity {proximity}
    zero_outside_proximity

    distortion_error {dist_error}
    default_point_error {point_error}

    coordinate_precision 5
    value_precision 4
    fill_grid
    calculate_control_point_stdres {calcstdres}
    columns lon lat de dn
    ''';

    df=result_file
    csr='true' if calcstdres else 'false'
    cmd=dgcmd
    cmd=cmd.replace('{data_file}',data_file)
    cmd=cmd.replace('{grid_spacing}','{0} {1}'.format(grid_spacing*1.5,grid_spacing))
    cmd=cmd.replace('{dist_error}',str(distortion_error))
    cmd=cmd.replace('{point_error}',str(point_error))
    cmd=cmd.replace('{proximity}',str(proximity))
    cmd=cmd.replace('{calcstdres}',csr)
    cmdf=df+'.cmd'
    with open(cmdf,'w') as dgf:
        dgf.write(cmd)
        if downgrades:
            for k,v in downgrades.iteritems():
                dgf.write("point {0} error {1}\n".format(k,v))
    call(('./datumgrid',cmdf,df))
    cpt=pd.read_csv(df+'_cpt.csv')
    grd=pd.read_csv(df+'_grd.csv')
    grdres=pd.read_csv(df+'_def.csv')
    return { 'cpt': cpt, 'grd': grd, 'grdres': grdres }

def write_log(text):
    with open(analysis_log,'a') as af:
        af.write(text)

def write_summary( version='',header=False):
    percentiles=(50.0,75.0,90.0,95.0,98.0,100.0)
    labels=["{0:.0f}".format(x) for x in percentiles]
    labels[-1]="max"
    if header or not version:
        with open(summary_file,'w') as f:
            f.write("version,")
            f.write(",".join('def'+x for x in labels))
            f.write(",")
            f.write(",".join('res'+x for x in labels))
            f.write(",err_gt20cm\n")
    if version:
        data=[version]
        dfm=pd.read_csv(calcdir+version+"_def.csv")
        data.extend("{0:.4f}".format(x) for x in np.percentile(dfm.distortion,percentiles))
        cpt=pd.read_csv(calcdir+version+"_cpt.csv")
        gt20=len(cpt[cpt.error > 0.2])
        data.extend(("{0:.4f}".format(x) for x in np.percentile(cpt.residual,percentiles)))
        data.append(str(gt20))
        with open(summary_file,'a') as f:
            f.write(",".join(data))
            f.write("\n")

print len(data)," points being used"
write_log("{0} control points being used\n\n".format(len(data)))

# Assess a suitable starting value for distortion error.  
# Decided that 5 was a good value as allowed most of the deformation to be absorbed...

if False:
    sumf=calcdir+'dg_summary.csv'
    with open(sumf,'w') as f:
        f.write("version,deferr,ssse,pc50,pc75,pc90,pc95,max\n")

    for i,error in enumerate((0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0,100.0)):
        error *= 5 # Changed point error from 0.01 to more realistic 0.05
        version='dg'+str(i)
        result=run_datumgrid(def_input_file,calcdir+version,distortion_error=error)
        cpt=result['cpt']
        ssse=np.sum(cpt.stdres*cpt.stdres)
        pc=np.percentile(cpt.residual,(50,75,90,95,100))
        with open(sumf,'a') as f:
            f.write("{7},{0},{1:.3f},{2:.3f},{3:.3f},{4:.3f},{5:.3f},{6:.3f}\n".format(
                error,ssse,pc[0],pc[1],pc[2],pc[3],pc[4],version))
    sys.exit()

# Look at rejecting nodes...
#
# Downweighting based on worst S.R.  
#  Iterative approach.  Downweight worst residual by factor of 10.  After doing this any other residuals within
#  a given radius are left till the next iteration.  Only top x% of S.Rs are considered.

spacing=np.exp(np.linspace(np.log(grid_spacing_initial),np.log(grid_spacing_final),max_iterations+1))
# Just need to sort out grid weighting so that finer grid doesn't increase distortion constraint...
# Experimentally it appears that reducing the spacing by a factor of 2 increases the influence of the
# distortion constraint, as there are 4 times as many grid points being summed.  So to counter this the 
# the distortion error must be increased by 2.
#
# 25.0 is intial grid spacing
weighting=base_distortion_error*(grid_spacing_initial/spacing)
weight_final=weighting[-1]

do_initial_reweighting = args.force_initial_reweighting
if not do_initial_reweighting and not args.ignore_initial_reweighting:
    do_initial_reweighting=not os.path.exists(calcdir+'final_grd.csv')

if do_initial_reweighting:
    print "Analysis phase 1: Initial station downweighting"
    downgrades={}

    
    result=run_datumgrid(
        def_input_file,
        calcdir+'reject0',
        grid_spacing=spacing[0],
        distortion_error=weighting[0],
        calcstdres=True
    )

    for i in range(max_iterations+1):
        write_log("\nInitial reweight iteration {0}\n".format(i))
        cpt=result['cpt']
        count=len(cpt.index)
        cpt.sort('stdres',ascending=False,inplace=True)
        cpt.index=range(count)
        nrej=int(count*percent_used/100.0)
        updated=False
        for nj in range(nrej):
            if cpt.stdres[nj] < min_std_res:
                break
            affected=False
            code=cpt.id[nj]
            for j in range(nj):
                dx=(cpt.lon[j]-cpt.lon[nj])*73000
                dy=(cpt.lat[j]-cpt.lat[nj])*100000
                ds=np.sqrt(dx*dx+dy*dy)
                if ds < influence_radius:
                    affected=True
                    break
            if affected:
                print "Skipping affected code:",code
                write_log("Skipping affected code: {0}\n".format(code))
                continue
            if code not in downgrades:
                downgrades[code]=0.01
            downgrades[code] /= downgrade_factor
            print code,"downgraded error to",downgrades[code]
            write_log("{0} downgraded error to {1:.5f}\n".format(code,downgrades[code]))
            updated=True
        # if not updated:
        #     break
        version='reject'+str(i)
        if i == max_iterations:
            version='reject_final'
        result=run_datumgrid(
            def_input_file,
            calcdir+version,
            downgrades=downgrades,
            grid_spacing=spacing[i],
            distortion_error=weighting[i],
            calcstdres=True
        )

    with open(calcdir+'downgrades.txt','w') as df:
        df.write("code,factor\n")
        for k in sorted(downgrades.keys()):
            df.write("{0},{1}\n".format(k,downgrades[k]))

# Now let's look at some iterative reweighting to get something more like a 1-norm.
# At each iteration we want to make sure that the total influence of the observations relative
# to the distortion is unchanged...
# Treat residuals less than min_residual as equivalent (ie don't think a fit better than this is 
# signficant - just good luck)
# Treat residuals greater than max_residual as disturbed, so downweight further

if True:
    print "Analysis phase 2: Iterative station reweighting for a range of distortion weightings" 

    version='reject_final'

    # Also want to reduce the distortion restriction somewhat to allow more locallised 
    # deformation, so do this in conjunction with reweighting.  Not sure if this really a good
    # idea!

    weights=(0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0)

    write_summary(header=True)
    if not args.ignore_initial_reweighting:
        write_summary('reject_final')
        finalgrid=pd.read_csv(calcdir+'reject_final_grd.csv')

    for w,weight in enumerate(weights):
        version='reweight'+str(w)
        write_log("\nIterative reweighting {0}: distortion error {1}\n".format(w,weight))
        # Start with result of final so bad stations are initially downweighted
        cptfile='' if args.ignore_initial_reweighting else 'final'
        lastgrid=None

        for i in range(reweight_iterations):
            write_log("  Iteration {0}\n".format(i))
            downgrades={}
            if cptfile:
                df=calcdir+cptfile+'_cpt.csv'
                cpt=pd.read_csv(df)
                err2=np.maximum(cpt.residual,min_residual)
                err2 *= np.maximum(1,cpt.residual/max_residual)
                for r in range(len(cpt.index)):
                    downgrades[cpt.id[r]]=err2[r]
            cptfile=version # Next iterations based on this file...
            result=run_datumgrid(
                def_input_file,
                calcdir+version,
                downgrades=downgrades,
                grid_spacing=grid_spacing_final,
                distortion_error=weight*weight_final
            )
            grid=result['grd']
            if lastgrid is not None:
                de=grid.de-lastgrid.de
                dn=grid.dn-lastgrid.dn
                ds=np.sqrt(np.max(de*de+dn*dn))
                write_log("    Maximum grid shift change {0:.4f}\n".format(ds))
                if ds < min_coord_change:
                    break
            lastgrid=grid

        write_summary(version)

        # Calc difference between grids

        if not args.ignore_initial_reweighting:
            rgrid=pd.read_csv(calcdir+version+'_grd.csv')
            dde=rgrid.de-finalgrid.de
            ddn=rgrid.dn-finalgrid.dn
            dds=np.sqrt(dde*dde+ddn*ddn)
            dif=pd.DataFrame({'lon':rgrid.lon,'lat':rgrid.lat,'dde':dde,'ddn':ddn,'dds':dds})
            dif.to_csv(calcdir+'final-'+version+'-grid.csv',index=False,cols=('lon','lat','dde','ddn','dds'))

# Recommend using version reweight4 as the residual errors match 
