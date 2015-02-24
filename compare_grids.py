import sys
import pandas as pd
import numpy as np
from subprocess import call

solutions=['final']
solutions.extend(['reweight'+str(i) for i in range(8)])

nicmodel='../../Deformation Model 20150101/model/patch_cc_2010_11/grid_cc_L1.csv'

percentiles=(50,75,90,95,98,100)

with open('../data/calcs/nic_comparison.csv','w') as cmpf:
    columns=['solution']
    columns.extend(('ds_pc'+str(p) for p in percentiles))
    columns[-1]='ds_max';
    columns.extend(('diff_pc'+str(p) for p in percentiles))
    columns[-1]='diff_max';
    cmpf.write(','.join(columns))
    cmpf.write("\n");
    for solution in solutions:
        finalmodel='../data/calcs/'+solution+'_grd.csv'
        diffmodel='../data/calcs/'+solution+'_nic_diff_grd.csv'

        call(('gridtool',
              'read','csv',finalmodel,
              'subtract','csv',nicmodel,
              'write','csv',diffmodel
             ))

        diff=pd.read_csv(diffmodel)
        diff['ds']=np.sqrt(diff.de*diff.de+diff.dn*diff.dn)

        fin=pd.read_csv(finalmodel)
        fin['ds']=np.sqrt(fin.de*fin.de+fin.dn*fin.dn)
        columns=[solution]
        columns.extend(('{0:.4f}'.format(x) for x in np.percentile(fin.ds,percentiles)))
        columns.extend(('{0:.4f}'.format(x) for x in np.percentile(diff.ds,percentiles)))
        cmpf.write(','.join(columns))
        cmpf.write("\n");

