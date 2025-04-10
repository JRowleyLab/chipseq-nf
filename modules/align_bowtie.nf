nextflow.enable.dsl=2

process ALIGN_BOWTIE2 {
    tag "$key"
    publishDir "${params.outdir}/alignment/$key/", pattern:'{*.log, *.bam}', mode: 'copy'

    input:
    tuple val(key), path(reads), val(input)
    path(index)

    output:
    tuple val(key), path("*bam"), val(input), emit: reads
    path("*.bowtie2.log"), emit: summary

    script:

        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/.rev.1.bt2//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/.rev.1.bt2l//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 \\
                -x \$INDEX \\
                --threads $params.threads \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                2> ${key}.bowtie2.log \\
                | samtools view -@ $params.threads -bhS -q 10 -o ${key}.bam -
                
        """  
}