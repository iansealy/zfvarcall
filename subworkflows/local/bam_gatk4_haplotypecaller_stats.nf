//
// Subworkflow to run GATK HaplotypeCaller and BCFtools stats
//

include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { BCFTOOLS_CONCAT       } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_INDEX        } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_STATS        } from '../../modules/nf-core/bcftools/stats/main'

workflow BAM_GATK4_HAPLOTYPECALLER_STATS {

    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_bai         // channel: [ val(meta), [ bai ] ]
    ch_fasta       // channel: [ val(meta), [ fasta ] ]
    ch_fasta_fai   // channel: [ val(meta), [ fai ] ]
    ch_fasta_dict  // channel: [ val(meta), [ dict ] ]
    ch_genome_bed  // channel: [ val(meta), [ bed ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: GATK HaplotypeCaller
    //
    ch_bam_bai_bed_1 = ch_bam.join(ch_bai, failOnDuplicate: true, failOnMismatch: true)
        .combine(ch_genome_bed)
        .map{ meta, bam, bai, meta2, bed -> [meta + meta2, bam, bai, bed, []] }
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_bed_1,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        [[], []], // no need for dbSNP file
        [[], []]  // no need for dbSNP index
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    //
    // MODULE: BCFtools concat
    //
    ch_vcfs_tbis = GATK4_HAPLOTYPECALLER.out.vcf.map{ meta, vcf -> [meta, vcf, meta.order] }
        .map { it[0].remove("order"); it }
        .groupTuple()
        .map{ meta, vcfs, order -> [meta, order.withIndex().sort().collect { vcfs[it[1]] }, []] }
    BCFTOOLS_CONCAT (
        ch_vcfs_tbis
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    //
    // MODULE: BCFtools index
    //
    BCFTOOLS_INDEX (
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    //
    // MODULE: BCFtools stats
    //
    ch_vcf_tbi = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_INDEX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
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
        BCFTOOLS_INDEX.out.tbi,
        BCFTOOLS_STATS.out.stats
    )

    emit:
    vcf      = BCFTOOLS_CONCAT.out.vcf         // channel: [ val(meta), [ vcf ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

// Dummy process

// Prevent nf-boost cleanup from prematurely deleting BAM files and indices
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
