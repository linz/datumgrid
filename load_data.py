#!/usr/bin/python

import sys
import fiona
from shapely.geometry import shape,Point
from shapely.prepared import prep
from csv import DictReader

# Load deformation model
modeldir='/home/ccrook/projects/deformation/models/published'
sys.path.append(modeldir+'/tools')
from LINZ.DeformationModel.Model import Model

components='patch_c1_20100904+patch_c2_20110222+patch_c3_20110613+patch_c4_20111223'

model=Model(modeldir+'/model',loadSubmodel=components)
model.setDate(2015.0,baseDate=2009.0)
ellipsoid=model.ellipsoid()

data_dir='../data/'
damage_file='LandDamage.shp'
point_file='eq_coords.csv'
model_data_file='model_data.csv'

damage_areas=[]
points=[]

with fiona.open(data_dir+damage_file) as df:
    for area in df:
        level=float(area['properties']['MaxEvent'])
        g=prep(shape(area['geometry']))
        damage_areas.append((level,g))

damage_areas.sort(key=lambda x:x[0], reverse=True)
#for d in damage_areas:
#    print 'Damage Level:',d[0]

with open(data_dir+model_data_file,'w') as mdf:
    mdf.write("code,lon,lat,de_obs,dn_obs,de_dm,dn_dm,de_dif,dn_dif,order,damage\n")
    with open(data_dir+point_file) as pf:
        dpf=DictReader(pf)
        for r in dpf:
            lon=float(r['lon_obs'])
            lat=float(r['lat_obs'])
            hgt=float(r['ellipsoid_hgt_obs'])
            rpt=Point(lon,lat)
            damage=0.0
            for level,shape in damage_areas:
                if shape.contains(rpt):
                    damage=level
                    break
            denu=model.calcDeformation(lon,lat)
            dedln,dndlt=ellipsoid.metres_per_degree(lon,lat,hgt)
            de_obs=(lon-float(r['lon_pre_eq']))*dedln
            dn_obs=(lat-float(r['lat_pre_eq']))*dndlt
            dif_de=de_obs-denu[0]
            dif_dn=dn_obs-denu[1]
            # Not all points have preeq heights
            # du_obs=hgt-float(r['ellipsoid_hgt_pre_eq'])
            order=max(r['coord_order_obs'],r['coord_order_pre_eq'])
            mdf.write("{0},{1:.7f},{2:.7f},{3:.4f},{4:.4f},{5:.4f},{6:.4f},{9:.4f},{10:.4f},{7},{8:.1f}\n".format(
                r['code'],lon,lat,de_obs,dn_obs,denu[0],denu[1],order,damage,dif_de,dif_dn))









