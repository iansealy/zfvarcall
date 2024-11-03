//
// Subworkflow to run freebayes, sort VCFs, index and generate stats
//

include { FREEBAYES       } from '../../modules/nf-core/freebayes/main'
include { BCFTOOLS_SORT   } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_STATS  } from '../../modules/nf-core/bcftools/stats/main'

workflow BAM_FREEBAYES_SORT_INDEX_STATS {

    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_bai         // channel: [ val(meta), [ bai ] ]
    ch_fasta       // channel: [ val(meta), [ fasta ] ]
    ch_fasta_fai   // channel: [ val(meta), [ fai ] ]
    ch_genome_bed  // channel: [ val(meta), [ bed ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: freebayes
    //
    ch_bam_bai_bed = ch_bam.join(ch_bai, failOnDuplicate: true, failOnMismatch: true)
        .combine(ch_genome_bed)
        .map{ meta, bam, bai, meta2, bed -> [meta + meta2, bam, bai, [], [], bed] }
    FREEBAYES (
        ch_bam_bai_bed,
        ch_fasta,
        ch_fasta_fai,
        [[], []], // no need for samples file
        [[], []], // no need for populations file
        [[], []]  // no need for copy number map file
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    //
    // MODULE: BCFtools sort
    //
    BCFTOOLS_SORT (
        FREEBAYES.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    //
    // MODULE: BCFtools concat
    //
    ch_vcfs_tbis = BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [meta, vcf, meta.order] }
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
        ch_bai,
        BCFTOOLS_CONCAT.out.vcf,
        BCFTOOLS_CONCAT.out.tbi,
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
    tuple val(meta2), path(bai)
    tuple val(meta3), path(vcf)
    tuple val(meta4), path(tbi)
    tuple val(meta5), path(stats)

    exec:
    Thread.sleep(1);
}
