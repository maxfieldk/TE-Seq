#use env: modbam
import pysam
import pandas as pd
import modbampy
from modbampy import ModBam
import matplotlib as mpl
import numpy as np
import yaml
import sys
import glob


print("script is running")
with open('conf/config.yaml', 'r') as file:
    conf = yaml.safe_load(file)

ref_annotation_dir = conf["reference_annotation_dir"]
rte_subfamily_read_level_analysis = conf["rte_subfamily_read_level_analysis"]
print("args")
print(sys.argv[0])
print(sys.argv[1])
print(sys.argv[2:])

reads_of_interest_path = sys.argv[1]
print("reads_of_interest_path", reads_of_interest_path)
bampaths = sys.argv[2:]
print("bampaths", bampaths)
samples = conf["samples"]
print("samples", samples)
sample_table = pd.read_csv(conf["sample_table"])
sample_table = sample_table.set_index('sample_name').loc[samples].reset_index()


element_df = pd.DataFrame(columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'type'])
for group in rte_subfamily_read_level_analysis:
    dftemp = pd.read_csv("%s/annotations/rte_beds/%s.bed"%(ref_annotation_dir, group), delimiter="\t")
    dftemp["type"] = group
    dftemp.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'type']
    element_df = pd.concat([element_df, dftemp])
element_dict = {}
for sample in samples:
    print(sample)
    element_dict[sample] = {}
    condition = sample_table.loc[sample_table["sample_name"] == sample, "condition"].iloc[0]
    bam_pattern = sample
    print(bam_pattern)
    pathtouse = "dog"
    for path in bampaths:
        print(path)
        if sample in path:
            print("match")
            pathtouse = path
    print("pathtouse is")
    print(pathtouse)
    for index, row in element_df.iterrows():
        print(row["chr"])
        reads = []
        with ModBam(pathtouse) as bam:
            for read in bam.reads(row[0], row[1], row[2]):
                for pos_mod in read.mod_sites:
                    ele = [*pos_mod]
                    reads.append(ele)
        try:
            df = pd.DataFrame(reads, columns=['id', 'refpos','querypos', 'modstrand', 'unclear','canon_base','mod_base', 'mod_base_score'])
            df["element_strand"] = row[5]
            df["element_start"] = row[1]
            df["element_stop"] = row[2]
            df["element_chr"] = row[0]
            df["element_type"] = row["type"]
            df["element_uid"] = "{}_{}_{}_{}".format(df["element_chr"][0], df["element_start"][0], df["element_stop"][0], df["element_strand"][0])
            df = df.assign(meth = np.where(df.mod_base_score > 255/2, 1, 0))
            df["condition"] = condition
            df["sample"] = sample
            element_dict[sample][df["element_uid"][0] + "_" + df["element_type"][0]] = df
        except:
            print("issue")

valuesextracted = [pd.concat(element_dict[sample].values(), axis=0) for sample in samples]
concatenated_dataframe = pd.concat(valuesextracted, axis=0)

concatenated_dataframe.to_csv(reads_of_interest_path, sep='\t', index=False)