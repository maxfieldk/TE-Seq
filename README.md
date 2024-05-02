# RTE Analysis

This project consists of a snakemake pipeline to analyze transposable element omics data. It consists of 4 modules, "Annotate Referene" (AREF), short-read RNA-Seq (SRNA), long-read RNA-Seq (LRNA), and long-read DNA-Seq (LDNA). LRNA and LDNA remain in active development, while AREF and SRNA are comparatively stable. Nevertheless changes to these modules will no doubt be forthcoming. A snakemake pipeline consists of several parts: a main "snakefile" which can deploy a number of module level snakefile, which contain rules. These rules take in inputs and produce outputs. In order to run this pipeline, a number of configuration files must be edited to reflect the end user's data and analytical desires. This README will walk you through the setup and runtime steps to run the AREF and SRNA modules. Accordingly, to get a particular "frozen" version of the pipeline, please specify a tag when cloning the pipeline (as shown below).
## Install snakemake
- Create a snakemake conda environment from which you can run the snakemake pipeline, and install the required snakemake executor plugins in your snamemake conda environment, e.g.
    ```
    mamba create --name snakemake snakemake snakemake-executor-plugin-slurm
    ```
- If you are not able to use docker/singularity and/or want to be able to modify environements, you can create the environments in the envs/ directory. E.g.:
   ```
   mamba env create --file envs/rseqc.yaml
   ```
   It will take time and occupy a substantial amount of disk space to recreate all of these environments. Using a container is the prefered way to deploy this pipeline.

## Setup your project directory
- Create a project directory
- Clone this pipeline into said directory, using a tag to specify a frozen version, or without one to get the latest version (this may give you more errors than a stable version).
  ```
  git clone -b v0.1.3 https://github.com/maxfieldk/RTE.git
  ```
- Copy the workflow/conf_example directory to ./conf
  ```
  cp workflow/conf_example conf
  ```
- Modify the contents of conf/sample_table_srna. Make sure sample names do not start with numbers; add an X in front if they do.
- Modify the contents of conf/config.yaml; make sure the samples, contrast, library_type are modified properly; also make sure that paths are present on your system and do not give you permissions errors. You will need to provide a number of standard genome annotation files, such as a gtf file, etc. This pipeline expects chromosome names to be in UCSC format ie. "chr1, chr2, ...".
  ```
  head {path}
  ```
- Create, in your project directory, the rawdata directory, and move your fastqs there. Make sure the naming is consistent with the naming scheme set forth in the conf/project_config_srna.yaml, which uses sample_name from the conf/sample_table_srna.csv i.e.:
```
source1: "rawdata/{sample_name}_R1.fastq.gz" source2: "rawdata/{sample_name}_R2.fastq.gz"
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

- peptable_srna.csv, is automatically updated each time you call snakemake, according to the rules set out in conf/project_config_srna.yaml applied to conf/sample_table.csv. This is how rawdata paths can be generated dynamically.
- Modify the contents of workflow/profile/default/config.yaml such that singularity containers have access to your data directory, or in general a directory which contains all files which are referenced in the pipeline and which contains your project directory.
  The workflow/profile/default/config.yaml instructs snakemake how to be run. More information on snakemake profiles can be found at https://snakemake.readthedocs.io/en/stable/executing/cli.html under the "Profiles" section header.
  ```
  singularity-args: '--bind /users/mkelsey/data'
  becomes
  singularity-args: '--bind /users/other_user/data'
  ```

- Performe a pipeline dry-run, which tells you which rules snakemake will deploy once really called
  ```
  conda activate snakemake
  #ensure you are in the pipeline directory which lives in your project folder, i.e. myproject/RTE/
  snakemake -n
  ```
- For help with snakemake, consult its highly usable and detailed docs at https://snakemake.readthedocs.io/en/stable/index.html
- For help with git, consult https://git-scm.com/docs/gittutorial

## Annotations
- T2T-hs1 annotations can be found at:
  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/current/
  Wherein *_T2T-CHM13v2.0_genomic.*.gz correspond to the refseq gtf or gff commonly used for RNA-Seq analysis.

