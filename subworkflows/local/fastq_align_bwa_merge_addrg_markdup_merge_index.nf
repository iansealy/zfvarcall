//
// Subworkflow to run BWA mem, Samtools merge (splits), GATK AddOrReplaceReadGroups,
// biobambam2 bamsormadup, Samtools merge (lanes) and Samtools index
//

include { BWA_MEM                                 } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SPLITS } from '../../modules/nf-core/samtools/merge/main'
include { GATK4_ADDORREPLACEREADGROUPS            } from '../../modules/nf-core/gatk4/addorreplacereadgroups/main'
include { BIOBAMBAM_BAMSORMADUP                   } from '../../modules/nf-core/biobambam/bamsormadup/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_LANES  } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                          } from '../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX {

    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_bwa_index  // channel: [ val(meta), [ bwa ] ]
    ch_fasta      // channel: [ val(meta), [ fasta ] ]
    ch_fasta_fai  // channel: [ val(meta), [ fai ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: BWA mem
    //
    BWA_MEM (
        ch_reads,
        ch_bwa_index,
        ch_fasta,
        true // Sort BAM file
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //
    // MODULE: Samtools merge (splits)
    //
    ch_bam_split = BWA_MEM.out.bam
        .map { it[0].remove("split"); it }
        .groupTuple()
    SAMTOOLS_MERGE_SPLITS (
        ch_bam_split,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SPLITS.out.versions)

    //
    // MODULE: GATK AddOrReplaceReadGroups
    //
    GATK4_ADDORREPLACEREADGROUPS (
        SAMTOOLS_MERGE_SPLITS.out.bam,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(GATK4_ADDORREPLACEREADGROUPS.out.versions)

    DUMMY1 (
        SAMTOOLS_MERGE_SPLITS.out.bam,
        GATK4_ADDORREPLACEREADGROUPS.out.bam
    )

    //
    // MODULE: biobambam2 bamsormadup
    //
    BIOBAMBAM_BAMSORMADUP (
        GATK4_ADDORREPLACEREADGROUPS.out.bam,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BIOBAMBAM_BAMSORMADUP.out.metrics.collect{it[1]})
    ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)

    DUMMY2 (
        GATK4_ADDORREPLACEREADGROUPS.out.bam,
        BIOBAMBAM_BAMSORMADUP.out.bam
    )

    //
    // MODULE: Samtools merge (lanes)
    //
    ch_bam_lanes = BIOBAMBAM_BAMSORMADUP.out.bam
        .map { it[0].remove("lane"); it[0].id = it[0].sample; it }
        .groupTuple()
    SAMTOOLS_MERGE_LANES (
        ch_bam_lanes,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_LANES.out.versions)

    //
    // MODULE: Samtools index
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_MERGE_LANES.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam      = SAMTOOLS_MERGE_LANES.out.bam    // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
    multiqc  = ch_multiqc_files                // channel: [ MultiQC files ]
}

// Dummy process

// Prevent nf-boost cleanup from prematurely deleting BAM files
process DUMMY1 {
    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bam)

    exec:
    Thread.sleep(1);
}
process DUMMY2 {
    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bam)

    exec:
    Thread.sleep(1);
}
