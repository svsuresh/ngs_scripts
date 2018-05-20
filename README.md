# RNA seq data analysis

RNA data can be used to call variants and/or for differential expression analysis. This repo contains bash scripts for running differential gene expression and variant calling. 
1. For variant calling, bash scripts for **GATK Best practices workflow for variant calling using RNAseq data** are provided. 
2. For differential gene expression analysis, **Tophat2-HTSEQ/Featurecounts** and **Tophat2-Cufflinks-Cummerbund (Tuxedo suite protocol 1)** scripts are provided. **cummerbund R script** is provided separate. 

Please note that tuxedo suite protocol 1 not suggested, rather new tuxedo workflow **HISAT2-StringTie-Ballgown**, is preferred in place of existing protocol. New protocol is provided in another repository of mine, as a [tuxedo 2 workflow](https://github.com/svsuresh/tuxedo2_snakemake) written using snakemake library. 
