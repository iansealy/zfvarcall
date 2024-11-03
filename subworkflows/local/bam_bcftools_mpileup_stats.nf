//
// Subworkflow to run BCFtools mpileup and generate stats
//

include { BCFTOOLS_MPILEUP  } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CONCAT   } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_STATS    } from '../../modules/nf-core/bcftools/stats/main'

workflow BAM_BCFTOOLS_MPILEUP_STATS {

    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_fasta       // channel: [ val(meta), [ fasta ] ]
    ch_genome_bed  // channel: [ val(meta), [ bed ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: BCFtools mpileup
    //
    ch_bam_bed = ch_bam.combine(ch_genome_bed)
        .map{ meta, bam, meta2, bed -> [meta + meta2, bam, bed] }
    BCFTOOLS_MPILEUP (
        ch_bam_bed,
        ch_fasta,
        false // no need to save mpileup
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

    //
    // MODULE: BCFtools concat
    //
    ch_vcfs_tbis = BCFTOOLS_MPILEUP.out.vcf.map{ meta, vcf -> [meta, vcf, meta.order] }
        .map { it[0].remove("order"); it }
        .groupTuple()
        .map{ meta, vcfs, order -> [meta, order.withIndex().sort().collect { vcfs[it[1]] }, []] }
    BCFTOOLS_CONCAT (
        ch_vcfs_tbis
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    //
    // MODULE: BCFtools stats
    //
    ch_vcf_tbi = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, vcf, tbi -> [meta, vcf, tbi] }
    BCFTOOLS_STATS (
        ch_vcf_tbi,
        [[], []], // no need for regions file
        [[], []], // no need for targets file
        [[], []], // no need for samples file
        [[], []], // no need for exons file
        [[], []]  // no need for fasta file
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    DUMMY (
        ch_bam,
        BCFTOOLS_MPILEUP.out.vcf,
        BCFTOOLS_MPILEUP.out.tbi,
        BCFTOOLS_STATS.out.stats
    )

    emit:
    vcf      = BCFTOOLS_CONCAT.out.vcf           // channel: [ val(meta), [ vcf ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

// Dummy process

// Prevent nf-boost cleanup from prematurely deleting BAM and VCF files
process DUMMY {
    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(vcf)
    tuple val(meta3), path(tbi)
    tuple val(meta4), path(stats)

    exec:
    Thread.sleep(1);
}
