import sys
import pandas as pd
import numpy as np
from scipy.stats import chi2


resfile, exposure, binary_outcome, outfile = sys.argv[1:5]

res = pd.read_csv(resfile, sep="\t")

if binary_outcome == "true":  # Different columns names for logistic regression
    res['logOR'] = np.log(res['OR'])
    res = res.rename(columns={'logOR': 'BETA', 'LOG(OR)_SE': 'SE'})

res = (res
       .rename(columns={'ID': 'SNPID', 'REF': 'Allele1', 'ALT':
                        'Allele2'})
       .drop_duplicates(subset=['SNPID', 'TEST'])
       .set_index(['SNPID', 'Allele1', 'Allele2'])
       .filter(['TEST', 'BETA', 'SE', 'P'])
       .query('TEST == "ADD" | TEST == "ADDx' + exposure + '"')
                .pivot(columns='TEST')
       .reset_index())
res.columns = ['SNPID', 'Allele1', 'Allele2',
               'Beta_Main', 'Beta_Interaction_1', 'SE_Beta_Main',
               'SE_Beta_Interaction_1_1', 'P_Value_Main', 'P_Value_Interaction']
res = (res
       .assign(Var_Beta_Main = lambda x: x.SE_Beta_Main ** 2,
               Var_Beta_Interaction_1_1 = lambda x: x.SE_Beta_Interaction_1_1
               ** 2)
       .assign(P_Value_Joint = lambda x: 1 - chi2.cdf(x.Beta_Main ** 2 /
                                                      x.Var_Beta_Main +
                                                      x.Beta_Interaction_1 ** 2
                                                      /
                                                      x.Var_Beta_Interaction_1_1,
                                                      df=2))
       .drop(['SE_Beta_Main', 'SE_Beta_Interaction_1_1'], axis='columns'))

res.to_csv(outfile, sep=" ", index=False, na_rep="NaN")
