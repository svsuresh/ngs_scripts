#For loop to declare each sample
for i in $(ls /home/suresh/rnaseq/raw_reads | grep ^[^d]| rev | cut -c 10- | rev | uniq)
do
# Make direcotry for fastqc
mkdir -p ${i%}fastqc
# Run fastqc
fastqc -o ./${i%}fastqc -f fastq /home/suresh/rnaseq/raw_reads/${i%}_R1.fastq /home/suresh/rnaseq/raw_reads/${i%}_R2.fastq
# Make directory to store cutadapt results
mkdir -p cutadapt
#Run cutadapt
cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC --quality-cutoff=20 --format=fastq -o ./cutadapt/${i%}_R1.cutadapt.fastq -p ./cutadapt/${i%}_R2.cutadapt.fastq /home/suresh/rnaseq/raw_reads/${i%}_R1.fastq /home/suresh/rnaseq/raw_reads/${i%}_R2.fastq
# Run top hat and save the out with directory extension _tophatout
tophat2 --no-coverage-search -o ./${i%}_tophatout -G /home/suresh/rnaseq/reference/genes_chr12.gtf -p 2 /home/suresh/rnaseq/reference/chr12 ./cutadapt/${i%}_R1.cutadapt.fastq ./cutadapt/${i%}_R2.cutadapt.fastq --rg-id ${i} --rg-sample ${i%} --rg-library rna-seq --rg-platform Illumina
# Index the bam
samtools index ./${i%}_tophatout/accepted_hits.bam
# Collect the stats
samtools flagstat ./${i%}_tophatout/accepted_hits.bam > ./${i%}_tophatout/accepted_hits.flastats
# Print only paired readsSort bam 
samtools view -bf 1 ./${i%}_tophatout/accepted_hits.bam -o ./${i%}_tophatout/accepted_hits.pe.bam 
# Run htseq count
htseq-count -f bam --stranded=no --mode=intersection-nonempty -r name ./${i%}_tophatout/accepted_hits.pe.bam /home/suresh/rnaseq/reference/genes_chr12.gtf > ./${i%}_tophatout/htseq.${i%}.txt
# Count features
featureCounts -B -p -t exon -g gene_id -a /home/suresh/rnaseq/reference/genes_chr12.gtf -o ./${i%}_tophatout/counts.txt ./${i%}_tophatout/accepted_hits.pe.bam
# make directory for storing htseq results
mkdir -p htseqresults
# Move all htseq results at one place
cp `ls ./${i%}_tophatout/htseq.${i%}.txt` htseqresults
done
# Run R script
Rscript /home/suresh/rnaseq/htseq.r
