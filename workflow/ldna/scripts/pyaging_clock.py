import pandas as pd
import pyaging as pya
from snakemake.script import snakemake

mydata = pd.read_csv(snakemake.input[0])
df = mydata.iloc[:,1:].T
adata = pya.pp.df_to_adata(df, imputer_strategy='knn')

pya.pred.predict_age(adata, ['Horvath2013', 'AltumAge', 'PCGrimAge', 'GrimAge2', 'DunedinPACE'])
adata.obs.to_csv(snakemake.output[0])