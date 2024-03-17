rule genomebrowserplots:
    input:
        rnasignalsF = expand("lrna/intermediates/{sample}/alignments/genome/{{alignmenttype}}/{sample}.F.bw", sample = config["lrna"]["samples"]),
        rnasignalsR = expand("lrna/intermediates/{sample}/alignments/genome/{{alignmenttype}}/{sample}.R.bw", sample = config["lrna"]["samples"]),
        dnamethylation = expand("ldna/intermediates/{sample}/methylation/{sample}_{{type}}_dss.tsv", sample = config["ldna"]["samples"])
    params:
        regions_of_interest = config["integrated"]["regions_of_interest"],
        outputdir = lambda w: "integrated/genomebrowserplots/%s"%w.alignmenttype
    output:
        "integrated/genomebrowserplots_{alignmenttype}.{type}.out"
    resources:
        cpus_per_task =10,
        mem_mb = 40000
    conda:
        "repeatanalysis"
    script: "scripts/genomebrowserplots.R" 

rule callgenomebrowserplots:
    input:
        expand("integrated/genomebrowserplots_{alignmenttype}.{type}.out", alignmenttype = ["dorado"], type = ["CG_m"])