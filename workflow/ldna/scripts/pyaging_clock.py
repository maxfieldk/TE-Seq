import pandas as pd
import pyaging as pya

mydata = pd.read_csv("/oscar/data/jsedivy/mkelsey/Nanopore/alz/RTE/cpg_formatted_methylation_filtered.csv")
df = mydata.iloc[:,1:].T
adata = pya.pp.df_to_adata(df, imputer_strategy='knn')

pya.pred.predict_age(adata, ['Horvath2013', 'AltumAge', 'PCGrimAge', 'GrimAge2', 'DunedinPACE'])
adata.obs.to_csv("pyaging_predictions.csv")