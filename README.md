# RTE-Seq: A Retrotransposable Element RNA-Seq Pipeline
This project consists of a __snakemake pipeline__ to analyze retrotransposable element (RTE) 'omics data.

To the unacquainted, the analysis of RTE, and more generally repetitive element, sequencing data can be a daunting task: the repetitive nature of these elements imposes  analytical pitfalls and raises a number of practical questions including:  
- Should I examine individual repetitive elements or rather groups of similar elements?
- Which elements and groups are most biologically significant?
- How do I deal with multi-mapping reads?

For these reasons, __repetitive elements are often neglected in RNA-sequencing analyses__.
This pipeline hopes to render investigation into these elements more tractable for those not steeped in the RTE literature.
It also aims address concerns pertaining to:
- non-referernce elements
- non-autonomous transcription driven by adjacent genes
- the quality of RTE annotations.

This project derives from my work in the __Sedivy Lab at Brown University__, where we study transposable elements, in particular __LINE1__.
#
![Asset 5](https://github.com/maxfieldk/RTE/assets/44215152/1f56451a-877b-4786-8303-3ead46a8a471)

## Pipeline Overview
  This pipeline conducts an end-to-end analysis of raw sequencing data, implementing state of the art RTE-minded computational methods. It produces a comprehensive analyses of repetitive element expression at both the level of an individual repetitive element as well as family groupings of these elements. It consists of 4 modules, "Annotate Referene" (AREF), short-read RNA-Seq (SRNA), long-read RNA-Seq (LRNA), and long-read DNA-Seq (LDNA). LRNA and LDNA remain in active development, while AREF and SRNA are comparatively stable.

  This pipeline will begin by fully annotating a user-provided genome for repetitive element content and identifying functionally important grouping of elements. This enriched annotation set allows the subsequent analysis to probe distinctions in expression between truncated, full-length, open reading-frame intact, RTEs as well as intergenic versus intragenic elements. Additionally, if provided with long-read DNA sequencing data derived from the specimens under investigation, this pipeline will call non-reference insertions using TLRD in order to create a custom genome and associated annotations describing all non-reference RTE element insertions. This enables the analysis of polymorphic RTE elements, which are likely to be the most active elements.  
    
  Then starting with raw sequencing reads (fastq files), this pipeline performs standard read level quality control, genomic alignment, and quantification of gene counts. Repeat-specific tools are then deployed to quantify repeat element expression by using both only uniquely mapped reads, or by using the information provided by unique mappings to guide probabilistic assignment of multi-mapping reads to specific repetitive elements. Differential transcript expression of repetitive and non-repetitive elements is assessed using DESeq2. Gene-set enrichment analyses are performed for both genes and families of repetitive elements. 
### Inputs
- Raw short-read RNA-Seq data stored in fastq file format.
- Optionally, long-read Nanopore DNA-Seq data
### Outputs
- Enhanced repetitive element annotations sets (gtf, bed, etc.)
- Normalized (and, if applicable, batch corrected) count tables
- Plots and figures
- If providing long-read Nanopore DNA sequences:
- A custom genome reference including contigs for non-reference RTE insertions
  Annotation sets which include these polymorphic RTEs  
- A phylogenetic and functional analysis of non-reference RTEs
### Limitations and various experimental considerations
- Before diving into implementation details and runtime instructions, some biological and technical considerations are in order. The assessment of RTE expression via short-read sequencing data has its limitations. Here I list several pitfalls to be aware of when designing experiments and interpreting results.
- Mappability is a function of read-length: a 50bp single-end read library will have a much greater number of multi-mapping reads as compared to a 150bp paired-end read library. Accordingly, insorfar as RTE expression analysis is a primary goal, I recommend opting for the latter.
- While overall a more accurate and less biased strategy than only considering uniquely mapping reads (cite), EM re-assignment of multi-mapping reads is imperfect, and element-level read estimates should be interpreted as best guesses when they are comprised of a substantial fraction of multi-mapping reads.
- Each individual human has only the order of 100 L1 insertions which are not captured by a reference genome. This speaks to the degree of polymorphism incurred by the RTE way of life. Analytically, this means that reads derived from such elements will necessarily be assigned to other, reference, elements. This problem is particularly acute seeing as the polymorphic, non-reference RTEs are likely some of the most active elements in the genome (they have had less time since insertion to accrue mutations versus older, and hence potentially fixed, elements). If nanopore DNA-sequences are not available, this problem can still be attenuated by using the most up to date reference genomes, such as telomere-to-telomere assemblies (cite), which provide a more accurate map of RTE insertions, especially in hard to assemble regions such as centromeres.
- If using an rRNA-depletion library preparation strategy, one should be aware that differences in RNA-quality (fragmentation)can have profound effects on repetitive element expression estimates. Degraded libraries will have a greater fraction of intronic reads, and seeing as introns contains many RTEs, our RTE signal will be artificially inflated. For this reason, poly-A selected libraries are preferred.
![Figure1_ReadDistribution_0714](https://github.com/maxfieldk/RTE/assets/44215152/ebca6a5a-1933-46f1-a505-d41965300d02)
- https://www.neb.com/en-us/products/e6310-nebnext-rrna-depletion-kit-human-mouse-rat
- The following twitter discussion may further illuminate the present concern: https://x.com/Faulkner_Lab/status/1578598533334458368
- Young L1 element sub-family (L1HS, L1PA2, etc.) transcripts contain a G-rich stretch of nucleotides in their 3' UTR which is thought to form a G-quadruplex secondary structure. Internal data show that during sequencing library amplification, can DNA polymerase can struggles with this structure, leading to artificially depressed count estimates. PCR free approaches such as long-read RNA sequencing may therefore provide a more accurate assessment of young L1 transcript abundance.
## Computational Requirements
  This pipeline uses the snakemake workflow manager, and consists of several parts: a main "snakefile" which can deploy a number of module level snakefiles, which in turn contain the rules which specify each step of the analysis. These rules are like functions, they take in inputs and produce outputs.  
### Software Requirements
  All software dependencies are packaged into a docker container which will automatically be built and used by the snakemake pipeline at the time of execution. Docker or singularity must be availible in order to have snakemake deploy containers. Users are free to alternatively choose to manually build the several conda environments required by the pipeline using the provided yaml environment specifications, though this alternative is likely to result in greater overall frustration. Snakemake itself is installed in a conda environment, and so users will need to have conda available regardless (link).  
  This pipeline was primarily developed on a compute cluster running a RedHat Linux OS and which uses the SLURM workload manager. Nevertheless the use of docker containers should enable users on other operating systems to run the pipeline without difficulty.  
### Hardware Requirements
  This pipeline allows for the parallel execution of many jobs which can occur simultaneously. Consequently, it is highly recommended to execute this pipeline on a compute cluster to take advantage of the corresponding diminishment of total runtime afforded by parallelization. Snakemake is designed to work with many commonly-used cluster workload managers such as SLURM.  
  Many steps require a substantial amount of RAM (north of 20 GB) to be available on your system, else they will fail and give you OOM (out-of-memory) errors.  
  The docker container is currently built around the X86 cpu architecture, and hence will not run on ARM based systems such as the latest apple silicon chips. A docker container suitable for these systems will be forthcoming in short order.  
# Installation
## Create a snakemake containing conda environment
  Create a snakemake conda environment from which you can run the snakemake pipeline, and install the required snakemake executor plugins in your snamemake conda environment, e.g.  
  ```
    conda create --name snakemake snakemake snakemake-executor-plugin-slurm
  ```
  If you are not able to use docker/singularity and/or want to be able to modify environements, you can create the environments in the envs/ directory. E.g.:  
  ```
   conda env create --file envs/rseqc.yaml
  ```
   It will take time and occupy a substantial amount of disk space to recreate all of these environments. Using a container is the prefered way to deploy this pipeline. I also recommend using the mamba command over the conda command if possible, environment building will proceed much more quickly
  ```
   mamba env create --file envs/rseqc.yaml
  ```
## Setup your project directory
  Create a project directory  
  Clone this pipeline into said directory, using a tag to specify a frozen version, or without one to get the latest version (this may give you more errors than a stable version).  
  ```
  git clone -b v0.1.6 https://github.com/maxfieldk/RTE.git
  ```
  Copy the workflow/conf_example directory to ./conf  
  ```
  cp workflow/conf_example conf
  ```
## Configure your analysis
  In order to run this pipeline, a number of configuration files must be edited to reflect your data and analytical decisions. Here I will walk you through the setup needed to execute the AREF and SRNA modules.
  
  Modify the contents of __conf/sample_table_srna.csv__. You provide two mandatory columns "sample_names" (e.g. Profiferating_1) and "condition" (e.g. Proliferating), and optionally meta-data variables, such as "batch" (which will be used to batch correct the differential expression analysis and provide batch-corrected counts, e.g. A). Make sure sample names do not start with numbers; add an X in front if they do.
  If you have added meta-data columns to your sample_table_srna which you would like to be represented as ordered factors in downstream visualizations (e.g. you are studying Alzheimer's samples and want Braak stage I to be followed by stages II, III, etc. in your plots), modify the contents of __conf/sample_table_source.R__ as suggested in the commented out code block of the script.

Modify the contents of __conf/config.yaml__. This file's contents determine the way in which the pipeline is run.

In the __"srna" section__, make sure that the values associated with the keys: "samples", "levels", "contrast", "library_type" are set. "samples" should be a list with every sample you plan to run the analysis on. Each value in "samples" and "levels" needs to have a corresponding entry in conf/sample_table_srna.csv, where "samples" corresponds to "sample_name", "levels" to "condition". "contrast" determines which contrasts are run in the differential expression analysis (e.g. if you have 9 samples spanning conditions A, B, and C, you can provide any or all of the following contrasts: "condition_B_vs_A", "condition_C_vs_A", "condition_C_vs_B". The condition which follows the _vs_ is the base comparator of the contrast, i.e. it is what the first condition is compared to. 

Then if you have genesets of interest, you can modify the "genesets_for_gsea" and "genesets_for_heatmaps" keys, adding a key-value pair for your geneset of interest. Please use gmt format for the former and a list of gene symbols in a text file for the latter (you can look at the provided genesets in the resources/ directory for further guidance).

In the __"aref" section__, locate the comment "#USER PROVIDED ANNOTATIONS - CHANGE PATHS AS NEEDED" and do as instructed. This set of annotations is all that you will need to provide to the pipeline, and consists of: (example urls to fetch the T2T-HS1 human reference are provided in parentheses)
- A reference genome (wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz)
- The repeatmasker.out file for this genome (wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz)
- Refseq gtf file (wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz)
- If using using nanopore reads to update your reference, a path to the repeatmasker executable on your system (https://www.repeatmasker.org/RepeatMasker/)
- Decompress all these files with the gunzip command, e.g.
  ```
  gunzip hs1.fa.gz
  ```
This pipeline expects chromosome names to be in UCSC format ie. "chr1, chr2, ...". You can easily inspect your annotation files to ensure they conform to this convention by using the following command:
  ```
  less {path_to_annotation}
  ```
  Annotations obtained from ncbi will typically need to have their naming convention changed.  
  To do so one need only run the chromToUcsc program, providing it with an annotation file and a chromAlias file. This is shown here: https://www.biostars.org/p/75369/#76420  
  ```
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chromAlias.txt
  chmod +x resources/programs/chromToUcsc
  chromToUcsc -g hs1 --get && chromToUcsc -i test.wig -o test.ucsc.wig -a hs1.chromAlias.tsv -g hs1
  ```
  I recommend you store these downloaded annotations in a directory one level above your project directory. Don't store any large files in resources/ seeing as it is under git control (git does not like large files). Just ensure that the paths are nested within a singularity bound directory if using the containerized workflow (read on to the segment on "workflow/profile/default/config.yaml" below for more details).  
  
Create, in your project directory, the srna/rawdata directory structure, and move your fastqs there. Make sure the naming is consistent with the naming scheme set forth in the conf/project_config_srna.yaml, which uses sample_name from the conf/sample_table_srna.csv i.e.:  
  ```
  source1: "srna/rawdata/{sample_name}_R1.fastq.gz" source2: "srna/rawdata/{sample_name}_R2.fastq.gz"
  ```
  If rawdata identifiers do not match with sample names, either you can rename them, or you can provide a mapping from sample_name to rawdata identifier. In this case you would add a column to your conf/sample_table_srna.csv titled something like  
  Then you would modify the "derive" block in your conf/project_config_srna.yaml as follows, chaning {sample_name} to your new column in conf/sample_table_csv, which in this case I titled {fq_sample_name}:  
  ```
  derive:
    attributes: [file_path_R1, file_path_R2]
    sources:
      source1: "srna/rawdata/{fq_sample_name}_R1_001.fastq.gz"
      source2: "srna/rawdata/{fq_sample_name}_R2_001.fastq.gz"
  ```
  peptable_srna.csv, is automatically updated each time you call snakemake, according to the rules set out in conf/project_config_srna.yaml applied to conf/sample_table.csv. This is how rawdata paths can be generated dynamically.  
  The workflow/profile/default/config.yaml instructs snakemake how to be run in your compute environment. More information on snakemake profiles can be found at https://snakemake.readthedocs.io/en/stable/executing/cli.html under the "Profiles" section header.  
  If using the containerized workflow, modify the contents of workflow/profile/default/config.yaml such that singularity containers have access to a directory which contains all files which are referenced in the pipeline and which contains your project directory.  
  ```
  singularity-args: '--bind /users/mkelsey/data,/oscar/data/jsedivy/mkelsey'
  becomes
  singularity-args: '--bind /users/YOURUSERNAME/data,/oscar/data/jsedivy/YOURUSERNAME'
  ```
## Workflow Logic:
### AREF
  In the spirit of use case flexibility, the AREF module has a number of workflow modifying parameters. These live in the aref section of the config.yaml file.  
  These parameters allow you to decide whether to create new annotations from scratch, to update annotations using long read sequencing data, or to use existing annotations which you had created during a previous run of the pipeline in another project (but which uses the same reference genome).  
  
  "symlink_aref" determines whether to run all the rules in aref module and thereby create an annotation set from scratch, or whether to merely symlink an existing directory of annotation files. If "symlink_aref" is set to "no", the "aref" module will be run (this is what you will do for your first analysis using this pipeline). If "symlink_aref" is set to "yes", the "aref" module will not be run, and the "aref" directory will be symlinked to the directory specified in "aref_dir" (this is recommended if you have previously completed an analysis with this pipeline, and with to re-use the annotations that analysis had created).  
  
  If you are not merely symlinking an existing AREF directory, you will need to specify whether you are creating an enhanced annotation set using long read DNA sequences (to call non-reference RTE insertions) or not.  
  The "update_ref_with_tldr" "response" value (yes or no) turn this feature on or off. If turning it on, you can specify whether to create one custom reference which all samples will use (e.g. if you had a number of long read sequencing data on spanning several conditions in ONE cell line) or to create one custom reference per sample (e.g. if you had long read sequencing done on multiple individuals). The "per_sample" key toggles between these two modes.  
  The "samples" , "sample_table", and "levels" keys in the aref section of the config are only relevant if creating a custom reference genome using long-read dna sequencing, and you can ignore these values if you are not using this feature. The same is true for the associated sample_table file, conf/sample_table_aref.csv, which will not be used in this case.  
### SRNA
  This module does not have workflow modifying parameters.  
## Deploying the pipeline
  First perform a pipeline dry-run this tells you which rules snakemake will deploy once really called  
  ```
  conda activate snakemake
  #ensure you are in the pipeline directory which lives in your project folder, i.e. myproject/RTE/
  snakemake -n
  ```
  If you are happy with this plan of action, deploy the pipeline by calling snakemake  
  ```
  snakemake
  ```
  I highly recommend familiarizing yourself with the basics of snakemake before embarking on a complex analysis with this pipelin. For help with snakemake, consult its highly usable and detailed docs at https://snakemake.readthedocs.io/en/stable/index.html  
  For help with git, consult https://git-scm.com/docs/gittutorial  
  If you encounter problems, please create a new issue on the github page.  
## Attribution
  The workflow graphic was adapted from work by Zandra Fagernas and carries a CC-BY 4.0 license.
