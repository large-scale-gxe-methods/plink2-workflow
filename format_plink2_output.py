import sys
import pandas as pd
from scipy.stats import chi2


resfile, exposure, outfile = sys.argv[1:4]

res = pd.read_csv(resfile, sep="\t")

res = (pd.read_csv(resfile, sep="\t")
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
