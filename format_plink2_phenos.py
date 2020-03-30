import sys
import pandas as pd

phenofile, sample_id_header, outcome, exposure, covar_names, delimiter, missing = sys.argv[1:8]

covars = [] if covar_names == "" else covar_names.split(" ")

output_cols = ["FID", sample_id_header, outcome, exposure] + covars
phenos = (pd.read_csv(phenofile, sep=delimiter, na_values=missing)
	.assign(FID=0)
	.loc[:, output_cols])
     
phenos.columns.values[:2] = ["FID", "IID"]
phenos.to_csv("plink2_phenotypes.txt", sep=" ", index=False, na_rep="NA")

num_covars = len(covars) + 1
plink_parameters_string = "1-" + str(num_covars + 2)
with open("plink2_parameters_string.txt", "w") as f:
	f.write(plink_parameters_string)
