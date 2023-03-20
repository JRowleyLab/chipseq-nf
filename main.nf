/* 
 * pipeline input parameters 
 */

//params.genome = "/Zulu/bnolan/Indexes/bwaIndex/hg38.fa"
//params.controlname = "control"
params.samplesheet = "${baseDir}/Samplesheets/samples_largetest_inputDependent2.csv"
params.outdir = "${baseDir}/results"
params.index = "/Zulu/bnolan/Indexes/Bowtie2Index/"
params.threads = "4"


log.info """\
         ===================================
         C H I P S E Q - N F   P I P E L I N E    
         ===================================
         outdir       : ${params.outdir}
         samplesheet  : ${params.samplesheet}
         threads      : ${params.threads}
         index        : ${params.index}

         """
         .stripIndent()

         

index_ch = Channel.value(file(params.index))


// Parse samplesheet and create reads channel            
Channel
        .from ( file(params.samplesheet) )
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ], row.input] }
        .set { read_pairs_ch }


// // Concatenate reads of same sample [e.g. additional lanes, resequencing of old samples]
//nfcore
read_pairs_ch
        .groupTuple()
        .branch {
            meta, fastq, input ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten(), input.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten(), input.flatten() ]
        }
        .set { ch_fastq }  


process CAT_FASTQ {
    //nf-core
    tag "$sample_id"
    publishDir "${params.outdir}/fastqMerged", mode: 'copy'

    input:
    tuple val(sample_id), path(reads, stageAs: "input*/*"), val(input) from ch_fastq.multiple

    output:
    tuple val(sample_id), path("*.merged.fastq.gz"), val(input) into cat_out_ch


    script:
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]

     if (readList.size >= 2) {
        def read1 = []
        def read2 = []
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }

        """
        cat ${read1.join(' ')} > ${sample_id}_1.merged.fastq.gz
        cat ${read2.join(' ')} > ${sample_id}_2.merged.fastq.gz
        """  
     }
}



//mix merged reads with singles
cat_out_ch
        .mix(ch_fastq.single)
        .into { cat_merged_ch; cat_merged_ch2 }



//  Run fastQC to check quality of reads files

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", pattern:"{*.html,fastqc_${sample_id}_logs}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(input) from cat_merged_ch

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}


// Trimming reads with Trim Galore

process trimming {
    tag "$key"

    publishDir "${params.outdir}/trimmed", pattern: "*.fq.gz", mode: 'copy'

    input: 
    tuple val(key), path(reads), val(input) from cat_merged_ch2

    output:
    tuple val(key), path("*.fq.gz"), val(input) into ch_out_trimmomatic, ch_out_trimmomatic2

    script:

    """
    trim_galore \\
                --paired \\
                ${reads[0]} \\
                ${reads[1]} \\
                --basename $key \\
                --cores 1 
    """
}

//  Run fastQC to check quality of reads files

process fastqc_trimmed {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", pattern:"{*.html,fastqc_${sample_id}_trimmed_logs}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(input) from ch_out_trimmomatic

    output:
    path("fastqc_${sample_id}_trimmed_logs") into fastqc_trimmed_ch

    script:
    """
    mkdir fastqc_${sample_id}_trimmed_logs
    fastqc -o fastqc_${sample_id}_trimmed_logs -f fastq -q ${reads}
    """  
}


// Align reads to index 

process align {
    tag "$key"
    publishDir "${params.outdir}/alignment/$key/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(reads), val(input) from ch_out_trimmomatic2
    path(index) from index_ch  

    output:
    tuple val(key), path("*bam"), val(input) into bam_ch, bam_ch2
    path("*.bowtie2.log"), optional: true into align_report_ch 

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


// Samtools stats for summary statistics on bwa-meth alignment
process samtools_stat_flagstat {
    tag "$key"
    publishDir "${params.outdir}/samtools/stats/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam), val(input) from bam_ch


    output:
    path('*stats') into samtools_stats_ch
    path('*flagstat') into samtools_flag_ch


    script:
    """
    samtools stats --threads ${params.threads} \\
                   $bam \\
                   > ${key}.stats

    samtools flagstat \\
                        --threads ${params.threads} \\
                        $bam \\
                        > ${key}.flagstat

    """
}



// samtools index
process samtools_sort {
    tag "$key"
    publishDir "${params.outdir}/samtools/bam/${key}/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam), val(input) from bam_ch2
    
    output:
    tuple val(key), path('*sort.bam'), val(input) into bam_sorted_channel

    script:
    """
    samtools sort $bam > ${key}.sort.bam
    """

} 


// Deduplicate bam files with picard
process deduplication{
    tag "$key"
    publishDir "${params.outdir}/picard/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam), val(input) from bam_sorted_channel 
    
    output:
    tuple val(key), path("*.deduplicated.bam"), val(input) into dedup_bam_ch, dedup_bam_ch2, dedup_bam_ch3
    path("*.MarkDuplicates.metrics.txt") into dedup_report_ch

    script:
    """
    picard \\
        MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${key}.deduplicated.bam \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=${key}.MarkDuplicates.metrics.txt
    """
}


//Combine replicates: split by '_', and group all samples
dedup_bam_ch
    .map { group_rep, bam, input ->
                        def(group) = group_rep.split("_")  
                        tuple( group, bam, input )
                        }
    .groupTuple()
    //Takes the first inputkey to avoid having nested tuples as inputkey
    .map { group, bam, input ->
                    tuple( group, bam, input.first() )

    }
    .into { dedup_groupsplit_ch; dedup_groupsplit_ch2 }



// Keep groups with more than 1 replicate, ready for combine replicates
dedup_groupsplit_ch
        .map {
            group, bams, input -> 
                        if (bams.size() != 1){ //Only keep 'groups' with >1 replicate
                            tuple( groupKey(group, bams.size()), bams, input)
                        }
        }
        .set{bam_sorted_groups_ch}



// Keep the individual replicates for inputs, if group has 1 replicate
dedup_groupsplit_ch2
        .map {
            group, bams, input -> 
                        if (bams.size() == 1 & group.contains('input')){ //Only keep 'groups' with >1 replicate
                            tuple( groupKey(group, bams.size()), bams, input)
                        }
                }
        .set{bam_single_inputs_ch}



// Combine replicates based on 'sample_rep' format, all 'sample' bam files will be merged
// Only combine replicates if there are replicates for that sample.
process combine_replicates {
    tag "$group"
    publishDir "${params.outdir}/samtools/merge/${group}/", pattern:'*', mode: 'copy'

    input:
    tuple val(group), path(bams), val(input) from bam_sorted_groups_ch

    output:
    tuple val(group), file("*.merged.bam"), val(input) into bam_merged_groups_ch

    script:
    """
    samtools merge -n ${group}.merged.bam $bams 
    """
}
  

//mix replicates with merged groups. If any inputs could be merged (per cell type), use the merged inputs for all further samples of that group
//Preferentially, only use merged inputs, if only 1 input is given for a cell type/ group, then use that one


dedup_bam_ch2 //deduplicated bam files
        .filter{ !it[0].contains('input') } //remove all input files
        .mix(bam_merged_groups_ch) //add merged inputs and IPs
        .mix(bam_single_inputs_ch) //add any individual replicates
        .set{bams_all_ch} //all input and IP bam files: tuple [sample, bam]


//CHECKPOINT

// samtools sort bam files
process samtools_sort2 {
    tag "$key"
    publishDir "${params.outdir}/samtools/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam), val(input) from bams_all_ch
    
    output:
    tuple val(key), path('*sort.bam'), val(input) into bams_sorted_ch3
    tuple val(key), path('*sort.bam') into bams_sorted_ch, bams_sorted_ch2

    
    script:
    """
    samtools sort $bam > ${key}.sort.bam
    """
}




process samtools_index {
    tag "$key"
    publishDir "${params.outdir}/samtools/bam/${key}/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bams_sorted_ch

    output:
    tuple val(key), path('*bai') into bam_indexed_ch

    script:

    """
    samtools index $bam
    """
}


// Add bai to bam channel for bamCoverage
bams_sorted_ch2
            .join(bam_indexed_ch)
            .set{bam_sorted_indexed_ch}
//CHECKPOINT




process bamCoverage {
    tag "$key"
    publishDir "${params.outdir}/bamCoverage/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam), path(bai) from bam_sorted_indexed_ch

    output:
    tuple val(key), path("*.bw") into bw_ch

    script:
    """
    bamCoverage -b $bam -o ${key}.bw --normalizeUsing BPM
    """
}

//CHECKPOINT
//MATCH INPUT TO EACH CONTROL USING THE 3RD ELEMENT IN EACH TUPLE [SAMPLE_ID, BAM, [INPUT]].
//NEW

bams_sorted_ch3 
        .branch{
            input: it[0].contains('input')
            ip: !it[0].contains('input')
        }
        .set {result}

//Match IP's to input's using the 3rd field.

result.ip
        .combine(result.input) 
        .map{
            ip, bam, ikey, input, bam2, na ->
                            def ipgroup = ip.minus(~/_.*/)
                            def inputkey = ikey.first()
                            def inputgroup = input.minus(~/_.*/)

                            if(inputkey.equals(inputgroup)){
                                tuple(ip, bam, bam2)
                                //tuple(ip, inputkey, input)
                            }
        }
        .set { ipbam_inputbam_ch }



// macs2

process macs2 {
    tag "$key"
    publishDir "${params.outdir}/macs2/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bamip), path(baminput) from ipbam_inputbam_ch

    output:
    tuple val(key), path("*.narrowPeak"), path(bamip) into peaks_bam_ch, peaks_bam_ch2, peaks_bam_ch3
    path("*.xls") into peaks_report_xls
    path("*")

    script:
    """
    macs2 \\
            callpeak \\
            -t $bamip \\
            -c $baminput \\
            -n $key \\
            -g hs \\
            --call-summits 
    """ 
}


peak_count_header_ch = Channel.fromPath("$projectDir/peak_count_header.txt", checkIfExists: true).toList()


process MULTIQC_CUSTOM_PEAKS {
    //nfcore/chipseq
    tag "$key" 

    input:
    tuple val(key), path(peak), path(bam) from peaks_bam_ch
    path(peak_count_header) from peak_count_header_ch

    output:
    path("*.peak_count_mqc.tsv") into custompeaks_count_multiqc_ch

    script:
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${key}", \$1 }' | cat $peak_count_header - > ${key}.peak_count_mqc.tsv
    """
}


process bamToBed {
    tag "$key" 

    input:
    tuple val(key), path(peak), path(bam) from peaks_bam_ch2

    output:
    tuple val(key), path(peak), path("*.bed") into bed_ch, bed_ch2 

    script:
    """
    bamToBed -i $bam > ${key}.bed
    """
}


// Create channel that joins 'control narrowPeak AND bam' with each treatment narrowPeak AND bam
// input 
// [T1, [T1.bam, T1.narrowPeak], [control.bam, control.narrowPeak]] ...


// Control
bed_ch
            .map {
                    group, peak, bed ->
                                if(group.contains('control')) {
                                    tuple( group, peak, bed )
                                } 
                    }
                    .set { control_peaks_bed_ch } 

// Treatments
bed_ch2
            .map {
                    group, peak, bed ->
                                if(!group.contains('control')) {
                                    tuple( group, peak, bed )
                                } 
                    }
                    .set { treatment_peaks_bed_ch } 

// Combine Treatment with Control
treatment_peaks_bed_ch
                    .combine(control_peaks_bed_ch)
                    .set {treat_cont_ch}

// Concept for allowing MAnorm to compare combinations of all samples. Need to handle this correctly to avoid huge numbers of comparisons.
// if (params.manormGroup){
//     // Control
//     bed_ch
//                 .map {
//                         group, peak, bed ->
//                                     if(group.contains('control')) {
//                                         tuple( group, peak, bed )
//                                     } 
//                         }
//                         .set { control_peaks_bed_ch } 

//     // Treatments
//     bed_ch2
//                 .map {
//                         group, peak, bed ->
//                                     if(!group.contains('control')) {
//                                         tuple( group, peak, bed )
//                                     } 
//                         }
//                         .set { treatment_peaks_bed_ch } 

//     // Combine Treatment with Control
//     treatment_peaks_bed_ch
//                         .combine(control_peaks_bed_ch)
//                         .set {treat_cont_ch}
// }else{
//     bed_ch
//         .combine(bed_ch2)
//         .map {
//             group, peak, bed, group2, peak2, bed2 ->
//                                             if(!group.contains(group2)) {
//                                                 tuple ( group,  group2 )
//                                             }
//         }
//         .set {treat_cont_ch}
// }

// manorm

process manorm{
    tag "$key_treatment"
    publishDir "${params.outdir}/manorm/${key_treatment}", pattern:"*", mode: 'copy'

    input:
    tuple val(key_treatment), path(peaks_treatment), path(bed_treatment), val(key_control), path(peaks_control), path(bed_control) from treat_cont_ch

    output:
    path('*_dir') into manorm_dir_ch

    script:
    """
    manorm \\
            --p1 $peaks_treatment \\
            --p2 $peaks_control \\
             --r1 $bed_treatment \\
             --r2 $bed_control \\
             --pf narrowpeak \\
             --rf bed \\
             -o ${key_treatment}v${key_control}_dir
    """ 
}


// Create multiqc report channel
  multiqc_ch = fastqc_ch
                    .mix(fastqc_trimmed_ch)
                    .mix(align_report_ch)                    
                    .mix(dedup_report_ch)
                    .mix(peaks_report_xls)
                    .mix(custompeaks_count_multiqc_ch)
                    .collect()


process multiqc {

    publishDir "${params.outdir}/multiqc/", pattern:"multiqc_report.html", mode:'copy'
       
    input:
    path('*') from multiqc_ch
    
    output:
    path('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 


workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong" )
}
