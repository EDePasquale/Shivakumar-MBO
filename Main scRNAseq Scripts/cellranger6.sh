module load cellranger/6.0.1

SAMPLE=${1}

cellranger count --id=${SAMPLE} --transcriptome=../../10x_reference/refdata-gex-GRCh38-2020-A --fastqs=../fastq/${SAMPLE} --localmem 128
