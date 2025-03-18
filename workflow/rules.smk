# rule genomebrowserplots:
#     input:
#         rnasignalsF = expand("lrna/intermediates/{sample}/alignments/genome/{{alignmenttype}}/{sample}.F.bw", sample = config["lrna"]["samples"]),
#         rnasignalsR = expand("lrna/intermediates/{sample}/alignments/genome/{{alignmenttype}}/{sample}.R.bw", sample = config["lrna"]["samples"]),
#         dnamethylation = expand("ldna/intermediates/{sample}/methylation/{sample}_{{type}}_dss.tsv", sample = config["ldna"]["samples"])
#     params:
#         regions_of_interest = config["integrated"]["regions_of_interest"],
#         outputdir = lambda w: "integrated/genomebrowserplots/%s"%w.alignmenttype
#     output:
#         "integrated/genomebrowserplots_{alignmenttype}.{type}.out"
#     resources:
#         cpus_per_task =10,
#         mem_mb = 40000
#     conda:
#         "repeatanalysis"
#     script: "scripts/genomebrowserplots.R" 

# rule callgenomebrowserplots:
#     input:
#         expand("integrated/genomebrowserplots_{alignmenttype}.{type}.out", alignmenttype = ["dorado"], type = ["CG_m"])


# rule epigenome_transcriptome_correlation:
#     input:
#         srna_results = "srna/results/agg/deseq/resultsdf.tsv",
#         lrna_results = "lrna/results/agg/deseq/dorado/relaxed/resultsdf.tsv",
#         ldna_methylation = expand("ldna/intermediates/{sample}/methylation/{sample}_CG_m_dss.tsv", sample = config["ldna"]["samples"]),
#         rteprommeth = "ldna/Rintermediates/perelementdf_promoters.tsv",
#         dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
#         dmls = "ldna/results/tables/dmls.CG_m.tsv"
#     output:
#         plots = "integrated/epigenome_transcriptome_correlation/objects/plots.rda"
#     resources:
#         cpus_per_task =10,
#         mem_mb = 60000
#     conda:
#         "ds"
#     script: "scripts/epigenome_transcriptome_correlation.R" 