##Condition for each sample
for i in $(ls ./raw_reads | grep ^[^d]| rev | cut -c 10- | rev | uniq)
do 

# Make direcotry for fastqc
#mkdir -p ${i%}_fastqc

# Run fastqc
#fastqc -o ./${i%}_fastqc -f fastq /home/suresh/rnaseq/raw_reads/${i%}_R1.fastq /home/suresh/rnaseq/raw_reads/${i%}_R2.fastq

# Make directory to store cutadapt results
#mkdir -p cutadapt

#Run cutadapt
#cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC --quality-cutoff=20 --format=fastq -o ./cutadapt/${i%}_R1.cutadapt.fastq -p ./cutadapt/${i%}_R2.cutadapt.fastq /home/suresh/rnaseq/raw_reads/${i%}_R1.fastq /home/suresh/rnaseq/raw_reads/${i%}_R2.fastq

# Run top hat and save the out with directory extension _tophatout
#tophat2 --no-update-check --no-coverage-search -o ./${i%}_tophat_out -G ./reference/genes_chr12.gtf -p 2 ./reference/chr12 ./cutadapt/${i%}_R1.cutadapt.fastq ./cutadapt/${i%}_R2.cutadapt.fastq --rg-id ${i} --rg-sample ${i%} --rg-library rna-seq --rg-platform Illumina

# Index the bam
#samtools index ./${i%}_tophat_out/accepted_hits.bam

# Collect the stats
#samtools flagstat ./${i%}_tophat_out/accepted_hits.bam > ./${i%}_tophat_out/accepted_hits.flastats

#Sort bam 
#samtools sort  -T /tmp/align.bam ./${i%}_tophat_out/accepted_hits.bam -o ./${i%}_tophat_out/accepted_hits.pos.sorted.bam  

# Cufflinks to aseemble transcripts and quantifies expressed features (genes, transcripts, TSS and CDS)
#cufflinks -g ./reference/genes_chr12.gtf ./${i%}_tophat_out/accepted_hits.pos.sorted.bam --no-update-check -b ./reference/chr12.fa -o ${i%}_cuff_out/ 

# Cuffcompare to generate comparable transcripts wrt reference transciptome
#mkdir -p MeOH_R3G_cuffcompare_out

cuffcompare -o ./MeOH_R3G_cuffcompare_out/cuffcompare.out -r ./reference/genes_chr12.gtf -s ./reference/chr12.fa ${i%}_cuff_out transcripts.gtf
find | grep transcripts.gtf$ > ./MeOH_R3G_cuffcompare_out/gtf.list

# Cuffmerge to generate master transcriptome file for all the samples including matching reference transcriptome 
mkdir -p .MeOH_R3G_cuffmerge_out/meoh_r3g
cuffmerge -o ./MeOH_R3G_cuffmerge_out/meoh_r3g -g ./MeOH_R3G_cuffcompare_out/cuffcompare.combined.gtf -s ./reference/chr12.fa ./MeOH_R3G_cuffcompare_out/gtf.list

# Cuffquant to quantify the counts
#cuffquant -o ./MeOH_R3G_cuffquant_out/ -p 2 --no-update-check -u -b ./reference/chr12.fa ./MeOH_R3G_cuffmerge_out/meoh_r3g/merged.gtf 

# Cuffdiff to output differential expression 
mkdir -p MeOH_R3G_cuffdiff_out
cuffdiff -L MeOH,R3G --multi-read-correct --no-update-check --FDR 0.1 -b ./reference/chr12.fa ./MeOH_R3G_cuffcompare_out/cuffcompare.combined.gtf -o ./MeOH_R3G_cuffdiff_out  ./MeOH_REP1_tophat_out/accepted_hits.pos.sorted.bam,./MeOH_REP2_tophat_out/accepted_hits.pos.sorted.bam,./MeOH_REP3_tophat_out/accepted_hits.pos.sorted.bam ./R3G_REP1_tophat_out/accepted_hits.pos.sorted.bam,./R3G_REP2_tophat_out/accepted_hits.pos.sorted.bam,./R3G_REP3_tophat_out/accepted_hits.pos.sorted.bam

done

# Cummerbund script
Rscript cummerbund.R
