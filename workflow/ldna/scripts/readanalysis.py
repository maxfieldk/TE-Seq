import logging
import pysam
import pandas as pd
import modbampy
from modbampy import ModBam
import matplotlib as mpl
import numpy as np
import yaml

print("script is running")
with open('conf/config.yaml', 'r') as file:
    conf = yaml.safe_load(file)

ref_annotation_dir = conf["reference_annotation_dir"]
rte_subfamily_read_level_analysis = conf["rte_subfamily_read_level_analysis"]

samples = conf["samples"]
sample_tablepath = conf["sample_table"]
sample_table = pd.read_csv(sample_tablepath)
sample_table = sample_table.set_index('sample_name').loc[samples].reset_index()

try:
    bampaths = snakemake.input["bams"]
    l1readspath = snakemake.output["l1reads"]
    samples = snakemake.params["samples"]
    print("snakemake variables loaded")
except:
    print("snakemake variables not loaded")
    samples = sample_table.sample_name
    l1readspath = 'results/tables/l1reads.tsv'
    bampaths = ["intermediates/%s/alignments/%s.sorted.filtered.bam"%(sample, sample) for sample in samples]

intactl1srefnonref = pd.read_csv(intactl1srefnonrefpath, delimiter="\t")

df = pd.DataFrame()
for group in rte_subfamily_read_level_analysis:
    dftemp = pd.read_csv("%s/annotations/rte_beds/%s.bed"%(ref_annotation_dir, group), delimiter="\t")
    dftemp.type = group
    df = pd.concat([df, dftemp])
    df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
##############

element_dict = {}
for sample in samples:
    element_dict[sample] = {}
    condition = sample_table.loc[sample_table["sample_name"] == sample, "condition"].iloc[0]
    bam_pattern = f"/*{sample_name}*/"
    bam_path = [path for path in inputs['bams'] if glob.fnmatch.fnmatch(path, pattern)][0]
    for index, row in df.iterrows():
        print(row["chr"])
        reads = []
        with ModBam(bam_path) as bam:
            for read in bam.reads(row["chr"], row["start"], row["end"]):
                for pos_mod in read.mod_sites:
                    ele = [*pos_mod]
                    reads.append(ele)
        try:
            df = pd.DataFrame(reads, columns=['id', 'refpos',
                                            'querypos', 
                                            'modstrand', 'unclear',
                                            'canon_base',
                                            'mod_base', 'mod_base_score'])
            df["element_strand"] = row["strand"]
            df["element_start"] = row["start"]
            df["element_stop"] = row["stop"]
            df["element_chr"] = row["chr"]
            df["element_type"] = row["family"]
            df["element_refstatus"] = row["source"]
            df["element_intactness"] = row["intactness"]
            df["element_uid"] = "{}_{}_{}_{}".format(df["element_chr"][0], df["element_start"][0], df["element_stop"][0], df["element_strand"][0])
            df = df.assign(meth = np.where(df.mod_base_score > 255/2, 1, 0))
            df["condition"] = condition
            df["sample"] = sample
            element_dict[sample][df["element_uid"][0] + "_" + df["element_type"][0]] = df
        except:
            print("exception", next(read.mod_sites))

valuesextracted = [pd.concat(element_dict[sample].values(), axis=0) for sample in samples]
concatenated_dataframe = pd.concat(valuesextracted, axis=0)

concatenated_dataframe.to_csv(l1readspath, sep='\t', index=False)