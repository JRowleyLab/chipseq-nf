# chipseq-nf
A Nextflow pipeline for the end-to-end data processing of ChIP-seq paired-end data. 

## How to use

1. Clone the repository: `git clone <repo.git>`
2. Create a singularity image from the nf-core Docker container: `singularity pull bsseq.sif docker://nolandocker/chipseq`
3. Create a samplesheet csv (`samples.csv`) that contains `sample`, `read 1` and `read 2` information like below:

**NOTE:** Use this exact format.

```
sample_id,fastq_1,fastq_2
Control_A,<path/to/read_1.fq.gz>,<path/to/read_2.fq.gz>
Control_B,<path/to/read_1.fq.gz>,<path/to/read_2.fq.gz>
```

4. Run the pipeline with the following code:

`nextflow run main.nf --samplesheet <samples.csv> --index <Bowtie2 index>`

## Replicates

The sample name and replicate must be separated by "_" (eg. **Sample_Replicate** OR **Control_A**). The pipeline will process single replicates AND combine replicates using `samtools merge`.

## Parameters

* `--samplesheet`: Samplesheet csv containing `sample_id,fastq_1,fastq_2`.
* `--outdir`: Output directory
* `--index`: path to index 

## Test

The pipeline can be tested with the following command `nextflow run main.nf --index <Bowtie2 index>`. Test reads are found in the `data/` folder.

## Singularity

The singularity image generated can be used with `-with-singularity <path/to/image>` or by placing the path in the `nextflow.config` file.
