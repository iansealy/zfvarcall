//
// Subworkflow to run fastp, FastQC and SeqKit split2
//

include { FASTP         } from '../../modules/nf-core/fastp/main'
include { FASTQC        } from '../../modules/nf-core/fastqc/main'
include { SEQKIT_SPLIT2 } from '../../modules/nf-core/seqkit/split2/main'

workflow FASTQ_FASTP_FASTQC_SPLIT {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: fastp
    //
    FASTP (
        ch_samplesheet,
        [],    // no need to specify adapters
        false, // no need to discard passing reads
        false, // no need to keep failed reads
        false  // no need to merge reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC (
        FASTP.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: SeqKit split2
    //
    SEQKIT_SPLIT2 (
        FASTP.out.reads
    )
    ch_reads_split = SEQKIT_SPLIT2.out.reads.transpose()
        .map{ meta, read -> [meta + [split: get_split_num(read.getName())], read]}
        .groupTuple()
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    DUMMY (
        FASTP.out.reads,
        SEQKIT_SPLIT2.out.reads
    )

    emit:
    reads    = ch_reads_split                  // channel: [ val(meta), [ reads ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
    multiqc  = ch_multiqc_files                // channel: [ MultiQC files ]
}

// Dummy process

// Prevent nf-boost cleanup from prematurely deleting reads
process DUMMY {
    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reads)

    exec:
    Thread.sleep(1);
}

// Helper function

// Extract split number from a read filename
// e.g. Get "2" from "ERS01_1.fastp.part_002.fastq.gz"
def get_split_num(String read) {
    return (read =~ /(\d+)\.(fastq|fq)(\.gz)?$/)[0][1].toInteger()
}
