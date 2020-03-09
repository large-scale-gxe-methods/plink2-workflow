import sys
import pandas as pd

phenofile, sample_id_header, outcome, covar_headers, exposure, delimiter, missing = sys.argv[1:8]
covars = covar_headers.split(" ")

covars.remove(exposure)

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