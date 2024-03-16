rule tldr:
    input:
        bam = lambda w: config["tldr_input_bam"]
    params:
        starting_ref = lambda w: config["starting_ref"]
    output:
        tldr = "tldr.table.txt"
    resources:
        cpus_per_task = 24,
        mem_mb = 100000,
        runtime = 300
    conda:
        "tldr"
    shell:
        """
mkdir -p tldr
cd tldr
tldr -b ../{input.bam} \
-e /oscar/data/jsedivy/mkelsey/tools/tldr/ref/teref.ont.human.fa \
-r {params.starting_ref} \
-p 20 \
--detail_output \
--extend_consensus 4000 \
--trdcol
cd ..
inputbam={input.bam}
mv tldr/${{inputbam%.*}}.table.txt {output.tldr}
        """

rule update_reference:
    input:
        reference = lambda w: config["starting_ref"],
        tldroutput = "updated_ref/tldr.table.txt"
    output:
        updated_reference = "updated_ref/updated_ref.fa"
    conda:
        "omics"
    script:
        "scripts/create_reference.R"

rule index_reference:
    input:
        reference = "updated_ref/updated_ref.fa"
    output:
        reference_index = "updated_ref/updated_ref.fa.fai"
    conda:
        "omics"
    shell:
        """
samtools faidx {input.reference}
        """
