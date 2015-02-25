import sys
import pandas as pd
import numpy as np
from subprocess import call

solutions=['final']
solutions.extend(['reweight'+str(i) for i in range(8)])

nicmodel='../../Deformation Model 20150101/model/patch_cc_2010_11/grid_cc_L1.csv'

percentiles=(50,75,90,95,98,100)

model_data={}
residual_data={}

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
        intmodel='../data/calcs/'+solution+'_nic_interp_grd.csv'

        call(('gridtool',
              'read','csv',finalmodel,
              'subtract','csv',nicmodel,
              'write','csv',diffmodel,
              'replace','csv',nicmodel,
              'write','csv',intmodel,
             ))

        diff=pd.read_csv(diffmodel)
        diff['ds']=np.sqrt(diff.de*diff.de+diff.dn*diff.dn)
        diff.to_csv(diffmodel,index=False)

        fin=pd.read_csv(finalmodel)
        fin['ds']=np.sqrt(fin.de*fin.de+fin.dn*fin.dn)
        columns=[solution]
        columns.extend(('{0:.4f}'.format(x) for x in np.percentile(fin.ds,percentiles)))
        columns.extend(('{0:.4f}'.format(x) for x in np.percentile(diff.ds,percentiles)))
        cmpf.write(','.join(columns))
        cmpf.write("\n");
        ds=[x for x in fin.ds]
        ds.sort();
        model_data[solution]=ds
        ds=[x for x in diff.ds]
        ds.sort();
        residual_data[solution]=ds

df=pd.DataFrame(data=model_data)
df.to_csv('../data/calcs/nic_model_shifts_sorted.csv',index=False)

df=pd.DataFrame(data=residual_data)
df.to_csv('../data/calcs/nic_differences_sorted.csv',index=False)


