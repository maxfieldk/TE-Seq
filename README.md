# TE-Seq: A Transposable Element Annotation & RNA-Seq Pipeline
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

This project consists of a __snakemake pipeline__ to analyze transposable element (TE) sequencing data.

To the unacquainted, the analysis of TE, and more generally repetitive element, sequencing data can be a daunting task: the repetitive nature of these elements imposes  analytical pitfalls and raises a number of practical questions including:  
- Should I examine individual repetitive elements or rather groups of similar elements?
- Which elements and groups are most biologically significant?
- How do I deal with multi-mapping reads?

For these reasons, __repetitive elements are often neglected in RNA-sequencing analyses__.
This pipeline hopes to render investigation into these elements more tractable for those not steeped in the TE literature.
It also aims address concerns pertaining to:
- non-refererence elements
- non-autonomous transcription driven by adjacent genes
- the quality of TE annotations.

This project derives from my work in the __Sedivy Lab at Brown University__, where we study transposable elements, in particular __LINE1__.
#
![workflow_graphic](https://github.com/user-attachments/assets/46fe2120-b7c9-46bb-84e8-a5b830a74c1d)


## Pipeline Overview
  This pipeline conducts an end-to-end analysis of raw sequencing data, implementing state of the art TE-minded computational methods. It produces a comprehensive analyses of repetitive element expression at both the level of an individual repetitive element as well as family groupings of these elements. It consists of 4 modules, "Annotate Referene" (AREF), short-read RNA-Seq (SRNA), long-read RNA-Seq (LRNA), and long-read DNA-Seq (LDNA). LRNA and LDNA remain in active development and will be formally released at a future date. Accordingly, this guide pertains only to the AREF and SRNA modules. 

  This pipeline begins by fully annotating a user-provided genome for repetitive element content and identifying biologically important grouping of elements. This enriched annotation set allows downstream analyses to probe distinctions in RNA expression between truncated, full-length, and open reading-frame intact TEs, as well as between intergenic versus intragenic TEs. Additionally, if provided with long-read DNA sequencing data derived from the specimens under investigation, this pipeline will call non-reference insertions using the TLDR program in order to create a custom genome and associated annotations describing all non-reference TE element insertions. This enables the analysis of polymorphic TE elements, which are likely to comprise some of the most active elements.  
    
  Once genome annotation is complete, starting with raw sequencing reads (fastq files) this pipeline performs standard read level quality control, genomic alignment, and quantification of gene counts. Repeat-specific tools are then deployed to quantify repeat element expression by using either only uniquely mapped reads (unique mode), or uniquely mapped reads + multi-mapping reads (multi mode). Both modes are run by default. Differential transcript expression of repetitive and non-repetitive elements is assessed using DESeq2. Gene-set enrichment analyses are performed for both genes and families of repetitive elements. 
### Inputs
- Reference genome and several associated genome annotations (more on this below; pipeline works with any animal species - tested on human, mouse, and rabbit).
- Short-read RNA-Seq data stored in fastq file format (can be either single-end or paired-end)
- Optionally, long-read Nanopore DNA-Seq data
### Outputs
- Enhanced repetitive element annotations sets (gtf, bed, etc.)
- Normalized (and, if applicable, batch corrected) count tables
- Plots and figures
- If providing long-read Nanopore DNA sequences:
- A custom genome reference including contigs for non-reference TE insertions
  Annotation sets which include these polymorphic TEs  
- A phylogenetic and functional analysis of non-reference TEs
### Limitations and various experimental considerations
- Before diving into implementation details and runtime instructions, some biological and technical considerations are in order. The assessment of TE expression via short-read sequencing data has its limitations. Here I list several pitfalls to be aware of when designing experiments and interpreting results.
- Mappability is a function of read-length: a 50bp single-end read library will have a much greater number of multi-mapping reads as compared to a 150bp paired-end read library. Accordingly, insorfar as TE expression analysis is a primary goal, I recommend opting for the latter.
- While overall a more accurate and less biased strategy than only considering uniquely mapping reads, EM re-assignment of multi-mapping reads is imperfect, and element-level read estimates should be interpreted as best guesses when they are comprised of a substantial fraction of multi-mapping reads. For more on this, see: Teissandier, Aurélie, Nicolas Servant, Emmanuel Barillot, and Deborah Bourc’his. “Tools and Best Practices for Retrotransposon Analysis Using High-Throughput Sequencing Data.”  _Mobile DNA_  10, no. 1 (December 29, 2019): 52.  [https://doi.org/10.1186/s13100-019-0192-1](https://doi.org/10.1186/s13100-019-0192-1).
- If your sequencing data contain batch effects, one must be mindful that batch effect correction is an imperfect procedure which can in itself introduce a bias. This problem is particularly acute when biological conditions of interest are imbalanced with respect to batch. This pipeline accounts for batch in two places. The first is in DESeq2's linear model; the second is when producing batch corrected count tables using limma. Note that this pipeline is limited in correcting for up to two variables of interest (e.g. batch and sex), as limma can handle at most two variables. DESeq2's linear modelling better handles the batch effect and it's p-value estimates (for individual genes and repetitive elements) will be more reliable than any statistics produced downstream of limma "batch-corrected" counts (comparisons between repetitive element family total counts). Therefore, in the context of an unbalanced dataset, TE family comparison's based on adjusted counts need to be interpreted with a serving of salt, and in the context of non-batch adjusted results. For a more thorough discussion of the issue, consider: Nygaard, Vegard, Einar Andreas Rødland, and Eivind Hovig. “Methods That Remove Batch Effects While Retaining Group Differences May Lead to Exaggerated Confidence in Downstream Analyses.” Biostatistics 17, no. 1 (January 1, 2016): 29–39. https://doi.org/10.1093/biostatistics/kxv027.
 
- Each individual human has only the order of 100 L1 insertions which are not captured by a reference genome. This speaks to the degree of polymorphism incurred by the TE way of life. Analytically, this means that reads derived from such elements will necessarily be assigned to other, reference, elements. This problem is particularly acute seeing as the polymorphic, non-reference TEs are likely some of the most active elements in the genome (they have had less time since insertion to accrue mutations versus older, and hence potentially fixed, elements). If Nanopore DNA-sequences are not available, this problem can still be attenuated by using the most up to date reference genomes, such as telomere-to-telomere assemblies (Nurk, Sergey, Sergey Koren, Arang Rhie, Mikko Rautiainen, Andrey V. Bzikadze, Alla Mikheenko, Mitchell R. Vollger, et al. “The Complete Sequence of a Human Genome.” Science 376, no. 6588 (April 2022): 44–53. https://doi.org/10.1126/science.abj6987), which provide a more accurate map of TE insertions, especially in hard to assemble regions such as centromeres.
- If using an rRNA-depletion library preparation strategy, one should be aware that differences in RNA-quality (fragmentation) can have profound effects on repetitive element expression estimates. Degraded libraries will have a greater fraction of intronic reads, and seeing as introns contains many TEs, our TE signal will be artificially inflated. For this reason, poly-A selected libraries are preferred. The following figure from NEB shows how different library preparation protocols affect the fraction of reads derived from different genomic compartments:
![Figure1_ReadDistribution_0714](https://github.com/maxfieldk/TE-Seq/assets/44215152/ebca6a5a-1933-46f1-a505-d41965300d02)
https://www.neb.com/en-us/products/e6310-nebnext-rrna-depletion-kit-human-mouse-rat
- The following article further illuminates the present concern: Faulkner, Geoffrey J. “Elevated L1 Expression in Ataxia Telangiectasia Likely Explained by an RNA-Seq Batch Effect.”  _NEURON_  111, no. 5 (March 1, 2023): 610–11.  [https://doi.org/10.1016/j.neuron.2023.02.007](https://doi.org/10.1016/j.neuron.2023.02.007).
- Young L1 element subfamily (L1HS, L1PA2, etc.) transcripts contain a G-rich stretch of nucleotides in their 3' UTR which is thought to form a G-quadruplex secondary structure. Internal data show that during sequencing library amplification, DNA polymerase can struggle with this structure, leading to artificially depressed count estimates. PCR free approaches such as long-read RNA sequencing may therefore provide a more accurate assessment of young L1 transcript abundance.
## Computational Requirements
  This pipeline uses the Snakemake workflow manager, and consists of several parts: a main "snakefile" which can deploy a number of module level snakefiles, which in turn contain the rules which specify each step of the analysis. These rules are like functions, they take in inputs and produce outputs.  
### Software Requirements
  All software dependencies are packaged into a Docker container which will automatically be built and used by the Snakemake pipeline at the time of execution. Docker or singularity must be available in order to have Snakemake deploy containers. Users are free to alternatively choose to manually build the several Conda environments required by the pipeline using the provided yaml environment specifications, though this alternative is likely to result in greater overall frustration. Snakemake itself is installed in a Conda environment, and so users will need to have Conda available regardless https://conda.io/projects/conda/en/latest/user-guide/getting-started.html.  
  This pipeline was primarily developed on a compute cluster running a RedHat Linux OS and which uses the SLURM workload manager. Nevertheless the use of docker containers should enable users on other operating systems to run the pipeline without difficulty.  
### Hardware Requirements
  This pipeline allows for the parallel execution of many jobs which can occur simultaneously. Consequently, it is highly recommended to execute this pipeline on a compute cluster to take advantage of the corresponding diminishment of total runtime afforded by parallelization. Snakemake is designed to work with many commonly-used cluster workload managers such as SLURM.  
  Many steps require a substantial amount of RAM (north of 20 GB) to be available on your system, else they will fail and give you OOM (out-of-memory) errors.  This pipeline will likely fail if run on a personal computer without at least 64 GB of RAM.
# Installation
## Create a Snakemake containing Conda environment
  Create a Snakemake Conda environment from which you can run the Snakemake pipeline, and install the required Snakemake executor plugins in your Snamemake Conda environment.  
  ```
    conda create --name snakemake snakemake snakemake-executor-plugin-slurm
  ```
  If you are not able to use Docker/Singularity and/or want to be able to modify environments, you can create all of the environments in the envs/ directory.
  ```
   conda env create --file envs/rseqc.yaml
   conda env create --file ...
   ...
  ```
   It will take time and occupy a substantial amount of disk space to recreate all of these environments. Using a container is the preferred way to deploy this pipeline. I also recommend using the equivalent Mamba command over the Conda command if possible - environment building will proceed much more quickly.
  ```
   mamba env create --file envs/rseqc.yaml
  ```
## Setup your project directory
  Create a project directory
  ```
  mkdir myproject
  ```
  Clone this pipeline into said directory, using a tag to specify a frozen version, or without one to get the latest version (this may give you more errors than a stable version).  
  ```
  cd myproject
  git clone https://github.com/maxfieldk/TE-Seq.git
  ```
  Copy the workflow/conf_example directory to ./conf  
  ```
cd TE-Seq
cp -r workflow/conf_example conf
  ```
## Configure your analysis
  In order to run this pipeline, a number of configuration files must be edited to reflect your data and analytical decisions. Here I will walk you through the setup needed to execute the AREF and SRNA modules.
  
  Modify the contents of __conf/sample_table_srna.csv__. You provide two mandatory columns "sample_names" (e.g. Profiferating_1) and "condition" (e.g. Proliferating), and optionally meta-data variables. You can include up to 2 columns which contain the string "batch" (e.g. "batch", and "batch_lane"), which will be used to batch correct the differential expression analysis and provide batch-corrected counts. If you only have one batch variable to model, it must be called "batch", and if you have two, the main batch effect should be titled "batch", and this will be the batch effect shown in various qc plots (but both will be modeled by Deseq2/limma). Make sure sample names do not start with numbers (add an X in front if they do), and that these names do not contain a period "." or dash "-" character (these are excluded by wildcard constraints in several rules).
  If you have added meta-data columns to your sample_table_srna which you would like to be represented as ordered factors in downstream visualizations (e.g. you are studying Alzheimer's samples and want Braak stage I to be followed by stages II, III, etc. in your plots), modify the contents of __conf/sample_table_source.R__ as suggested in the commented out code block of the script.

Modify the contents of __conf/config.yaml__. This file's contents determine the way in which the pipeline is run.

In the __"srna" section__, make sure that the values associated with the keys: "samples", "levels", "contrast", "library_type" are set. "samples" should be a list with every sample you plan to run the analysis on. Each value in "samples" and "levels" needs to have a corresponding entry in conf/sample_table_srna.csv, where "samples" corresponds to "sample_name", "levels" to "condition". "contrast" determines which contrasts are run in the differential expression analysis (e.g. if you have 9 samples spanning conditions A, B, and C, you can provide any or all of the following contrasts: "condition_B_vs_A", "condition_C_vs_A", "condition_C_vs_B". The condition which follows the _vs_ is the base comparator of the contrast, i.e. it is what the first condition is compared to (e.g. condition_treatment_vs_control). 

Then if you have genesets of interest, you can modify the "genesets_for_gsea" and "genesets_for_heatmaps" keys, adding a key-value pair for your geneset of interest. Please use gmt format for the former and a list of gene symbols in a text file for the latter (you can look at the provided genesets in the resources/ directory to see exactly how these files are to be formatted).

In the __"aref" section__, locate the comment "#USER PROVIDED ANNOTATIONS - CHANGE PATHS AS NEEDED". This set of annotations is all that you will need to provide to the pipeline, and consists of several files. You will have to obtain these files, and specify their paths in the config.yaml aref section. Here I show how to do this for the T2T-HS1 human genome or the MM39 mouse genome, but analogous files can be found for most animal species.
***
***
### HUMAN - Genome: T2T-HS1
- Download these files, and place them in a genome directory adjacent to your project directory
```
cd ..
mkdir genome_files
cd genome_files
#reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
#repeatmasker.out file (technically optional if you want the pipeline to repeatmask for you: change the "run_repeatmasker" "response" to "yes"; but beware! this is a computationally expensive step)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz
#Refseq gtf file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
#Refseq gff3 file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
```
- These files will need to be decompressed
```
gunzip hs1.fa.gz; mv hs1.fa reference.ucsc.fa
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz; mv GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf refseq.gtf
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz; mv GCF_009914755.1_T2T-CHM13v2.0_genomic.gff refseq.gff3
gunzip hs1.repeatMasker.out.gz; mv hs1.repeatMasker.out repeatmasker.ucsc.out
```
This pipeline expects chromosome names to be in UCSC format ie. "chr1, chr2, ...". You can easily inspect your annotation files to ensure they conform to this convention by using the following command:
  ```
  less {path_to_annotation}
  ```
  Annotations obtained from NCBI will typically need to have their naming convention changed.  
  To do so one need only run the resources/programs/chromToUcsc program, providing it with a sorted annotation (or genome reference) file and a chromAlias file. This is shown here: https://www.biostars.org/p/75369/#76420  
  ```
  sort -k1,1V -k4,4n -k5,5n refseq.gtf > refseq.sorted.gtf
  sort -k1,1V -k4,4n -k5,5n refseq.gff3 > refseq.sorted.gff3
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chromAlias.txt
  chmod +x ../TE-Seq/resources/programs/chromToUcsc
  ../TE-Seq/resources/programs/chromToUcsc -i refseq.sorted.gtf -a hs1.chromAlias.txt > refseq.sorted.ucsc.gtf
  ../TE-Seq/resources/programs/chromToUcsc -i refseq.sorted.gff3 -a hs1.chromAlias.txt > refseq.sorted.ucsc.gff3
rm refseq.gtf; rm refseq.sorted.gtf; rm refseq.gff3; rm refseq.sorted.gff3
  ```
***
***
### MOUSE - Genome: MM39
- Download these files, and place them in a genome directory adjacent to your project directory
```
cd ..
mkdir genome_files
cd genome_files
#reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
#repeatmasker.out file (technically optional if you want the pipeline to repeatmask for you: change the "run_repeatmasker" "response" to "yes"; but beware! this is a computationally expensive step)
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.out.gz
#Refseq gtf file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
#Refseq gff3 file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz
```
- These files will need to be decompressed
```
gunzip mm39.fa.gz; mv mm39.fa reference.ucsc.fa
gunzip GCF_000001635.27_GRCm39_genomic.gtf.gz; mv GCF_000001635.27_GRCm39_genomic.gtf refseq.gtf
gunzip GCF_000001635.27_GRCm39_genomic.gff.gz; mv GCF_000001635.27_GRCm39_genomic.gff refseq.gff3
gunzip mm39.out.gz; mv mm39.out repeatmasker.ucsc.out
```
This pipeline expects chromosome names to be in UCSC format ie. "chr1, chr2, ...". You can easily inspect your annotation files to ensure they conform to this convention by using the following command:
  ```
  less {path_to_annotation}
  ```
  Annotations obtained from ncbi will typically need to have their naming convention changed.  
  To do so one need only run the resources/programs/chromToUcsc program, providing it with a sorted annotation (or genome reference) file and a chromAlias file. This is shown here: https://www.biostars.org/p/75369/#76420  
  ```
  sort -k1,1V -k4,4n -k5,5n refseq.gtf > refseq.sorted.gtf
  sort -k1,1V -k4,4n -k5,5n refseq.gff3 > refseq.sorted.gff3
  wget http://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/chromAlias.txt.gz
  gunzip chromAlias.txt.gz
  chmod +x ../TE-Seq/resources/programs/chromToUcsc
  ../TE-Seq/resources/programs/chromToUcsc -i refseq.sorted.gtf -a chromAlias.txt > refseq.sorted.ucsc.gtf
  ../TE-Seq/resources/programs/chromToUcsc -i refseq.sorted.gff3 -a chromAlias.txt > refseq.sorted.ucsc.gff3
  rm refseq.gtf; rm refseq.sorted.gtf; rm refseq.gff3; rm refseq.sorted.gff3
  ```
***
***
  I recommend you store these downloaded annotations in a directory one level above your project directory (place genomes_files/ adjacent to the TE-Seq/ directory). Just ensure that the paths are nested within a Singularity bound directory if using the containerized workflow (read on to the segment on "workflow/profile/default/config.yaml" below for more details).  

Create, in your project directory, the srna/rawdata directories, and move your fastqs there. 
```
mkdir -p srna/rawdata
#now move all your fastqs into this directory
```
Make sure your fastq file naming is consistent with the naming scheme set forth in the conf/project_config_srna.yaml. This scheme fills in the sample_name from the conf/sample_table_srna.csv for both the read1 and read2 fasts file paths i.e.:  
  ```
  source1: "srna/rawdata/{sample_name}_R1.fastq.gz" source2: "srna/rawdata/{sample_name}_R2.fastq.gz"
  ```
  If rawdata file names do not match with sample names, either you can rename them, or you can provide a mapping from sample_name to rawdata identifier. In this case you would add a column to your conf/sample_table_srna.csv titled something like  "fq_sample_name", and add the rawdata identifiers to this column. Then you would modify the "derive" block in your conf/project_config_srna.yaml as follows, changing {sample_name} to your new column name in conf/sample_table_csv, as done below "fq_sample_name":  
  ```
  derive:
    attributes: [file_path_R1, file_path_R2]
    sources:
      source1: "srna/rawdata/{fq_sample_name}_R1_001.fastq.gz"
      source2: "srna/rawdata/{fq_sample_name}_R2_001.fastq.gz"
  ```
  peptable_srna.csv, is automatically updated each time you call Snakemake, according to the rules set out in conf/project_config_srna.yaml applied to conf/sample_table.csv. This is how rawdata paths can be generated dynamically.  
  The workflow/profile/default/config.yaml instructs Snakemake how to be run in your compute environment. More information on Snakemake profiles can be found at https://snakemake.readthedocs.io/en/stable/executing/cli.html under the "Profiles" section header.  
  If using the containerized workflow (default), modify the contents of workflow/profile/default/config.yaml such that Singularity containers have access to a directory which contains all files which are referenced in the pipeline and which contains your project directory.  
  ```
  SUPPOSE YOUR PIPELINE DIRECTORY'S PATH IS /users/YOURUSERNAME/data/myproject/TE-Seq
  AND THAT YOU HAVE SUPPLIED GENOME ANNOTATIONS WITH THE FOLLOWING PATH /oscar/data/jsedivy/YOURUSERNAME/refseq.gtf
  CHANGE
  singularity-args: '--bind /users/mkelsey/data,/oscar/data/jsedivy/mkelsey'
  TO
  singularity-args: '--bind /users/YOURUSERNAME/data,/oscar/data/jsedivy/YOURUSERNAME'
  ```
## Workflow Logic:
### AREF
  In the spirit of use case flexibility, the AREF module has a number of workflow modifying parameters. These live in the aref section of the config.yaml file.  
  These parameters allow you to decide whether to create new annotations from scratch, to update annotations using long read sequencing data, or to use existing annotations which you had created during a previous run of the pipeline in another project (but which uses the same reference genome).  
  
  "symlink_aref" determines whether to symlink (similar to copying, but without duplicating the files) an existing directory of annotation files, or if instead all the rules in the aref module should be run thereby creating an annotation set from scratch. If "symlink_aref" "response" is set to "no", the "aref" module will be run (this is the default). If "symlink_aref" is set to "yes", the "aref" module will not be run, and the "symlink_aref" "aref_symlinkdir" will be symlinked to the "aref" (this is recommended if you have previously completed an analysis with this pipeline, and with to re-use the annotations that analysis had created).   
  
  If you are not symlinking an existing AREF directory, you will need to specify whether you are creating an updated genome by providing long-read Nanopore DNA sequences (to call non-reference TE-Seq insertions) or not.  
  The "update_ref_with_tldr" "response" value (yes or no) turns this feature on or off. If turning it on, you can specify whether to create one custom reference which all samples will use (e.g. if you had a number of long read sequencing data on spanning several conditions in ONE cell line) or to create one custom reference per sample (e.g. if you had long read sequencing done on multiple individuals). The "per_sample" key toggles between these two modes.  
  The "samples" , "sample_table", and "levels" keys in the aref section of the config are only relevant if creating a custom reference genome using long-read dna sequencing, and you can ignore these values if you are not using this feature. The same is true for the associated sample_table file, conf/sample_table_aref.csv, which will not be used in this case. 
### SRNA
  The "per_sample_ref" key's value instructs the pipeline as to whether each sample has its own unique reference and TE annotations ("yes"; this situtaion occurs when you have nanopore DNA sequencing on all samples and wish to use custom references), or whether all samples will be using the same reference and TE annotations.  
  In order to rule out mycoplasma contamination of cultured cells, reads are mapped to a collection of mycoplasma genomes (nearly 100% of reads should then fail to map). This is done by default. Insofar as your samples are not cultured cells, you can omit this step by switching the "srna" "map_to_mycoplasma" key's value from "yes" to "no". 
## Applying this pipeline to non-human/mouse species.
Several additional parameters will need to be specified.
1. Change the aref species key's value in conf/config.yaml to a valid NCBI Taxonomy Database species name which is also contained in the RepeatMasker repeat database. 
2. If a repeatmasker.out file is not availible for your species, you can have the pipeline run repeatmasker and generate this file for you. Merely change the value of the "run_repeatmasker" "response" to "yes" (at which point the starting_ref_repeatmasker path is ignored). Note that genome-wide repeat masking is a computationally expensive step to be avoided if possible.
3. In the event that your gene annotation gtf's gene_id column does not contain symbol identifiers, you can provide a mapping from whatever identifier the gtf currently hs a(possibly an ENTREZ ID number) to a gene symbol (e.g. TP53). This is essential for the enrichment analysis by GSEA to work. Change the "gtf_id_mapping" "response" to "yes" and provide a "path" to a tsv file with two columns named "gene_id" and "gene_name" - gene_id will be swapped for gene_name in during enrichment analysis. 
4. If using long DNA reads to call novel insertions, specify a fast file in the "aref" "tldr_te_ref" section of conf/config.yaml . This key's value should be a file path to a fasta file containing consensus sequences for all TEs which will be probed for novel insertions by the TLDR program. The naming convention for each fasta element should be akin to that found in the mouse/human files found in the resources directory, i.e. ">famlily:element", or for a concrete example: ">L1:L1HS". The family and element names should be ones found in your repeatmasker annotation.
Optionally
5. In the aref/scripts/annotate_rtes.R script:
CTRL F to the only occurrence of "conf$species". You will find two blocks of code which provide species specific annotations for repeat families which are extensively examined in downstream scripts. In the same manner as is done in these two blocks for human and mouse: add an else {...} block, and supply your own regex for whatever repeat groupings you are interested in. If you do not do this, the script will automatically select the least diverged LTR, SINE, and LINE families for downstream scrutiny.
## Deploying the pipeline
  First perform a pipeline dry-run. This tells you which rules Snakemake will deploy once called "for real".   Ensure you are calling this command from within your pipeline directory (TE-Seq) which lives in your project folder, i.e. myproject/TE-Seq/
  ```
  conda activate snakemake
  snakemake --profile workflow/profile/default -n
  ```
  If you are happy with this plan of action, deploy the pipeline by calling Snakemake  
  ```
  snakemake --profile workflow/profile/default
  ```
  If all rules are executed successfully, you will not recieve an error message. At this point you can run the following command to generate a report which consolidates some of the most insightful visualizations. Note that not every plot is included in this report, seeing as the resulting html file would be too large for practical use.
  ```
  snakemake --profile workflow/profile/default --report report.html
  ```
  In the event that a rule fails, all hope is not lost!
  I recommend familiarizing yourself with the basics of Snakemake before embarking on a complex analysis with this pipeline. For help with Snakemake, consult its highly usable and detailed docs at https://snakemake.readthedocs.io/en/stable/index.html  
  For help with git, consult https://git-scm.com/docs/gittutorial  
  If you encounter problems, please create a new issue on the github page.  


## Common Issues

```
r$> conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
Error in configr::read.config(file = "conf/config.yaml")[[module_name]] : 
  subscript out of bounds
  ```
This error likely indicates that your config.yaml file is invalid (perhaps a dupliate key). You can copy the contents of your yaml to an online yaml checker to determine whether this is the cause, e.g. https://yamlchecker.com

---
```
ERROR: Failed to open file: "PATH"
```
An error of this type, assuming the path indeed exists and you have requisite permissions, could be due to a failure of having bound paths correctly in containerized rules. If the path you are specifying is to for instance, rawdata, ensure that the workflow/profile/default/config.yaml  "singularity-args" key specifies a bind path via which this path is accessible (i.e. has at least one path to a directory in which is nested the path of interest). If it does not, you can append a path to the directory which contains the path of interest. 

---
## Attribution
  Pending the publication of our manuscript, if you use this pipeline or any of its component parts, please cite this github page. 
    The workflow graphic contains icons which were adapted from work by Zandra Fagernas. It carries a CC-BY 4.0 license.  
