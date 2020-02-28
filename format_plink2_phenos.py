import sys
import pandas as pd

phenofile, sample_id_header, outcome, covar_headers, exposure, delimiter, missing = sys.argv[1:8]
covars = covar_headers.split(” “)

covars.remove(exposure)

output_cols = [sample_id_header, sample_id_header, outcome, exposure] + covars
phenos = (pd.read_csv(phenofile, sep=delimiter, na_values=missing)
     .loc[:, output_cols])
     
phenos.columns.values[:2] = [“FID”, “IID”]
phenos.to_csv(“plink2_phenotypes.csv”, sep=” “, index=False, na_rep=“NA”)