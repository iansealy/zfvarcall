//
// Subworkflow to run freebayes, sort VCFs, index and generate stats
//

include { FREEBAYES      } from '../../modules/nf-core/freebayes/main'
include { BCFTOOLS_SORT  } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX    } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'

workflow BAM_FREEBAYES_SORT_INDEX_STATS {

    take:
    ch_bam       // channel: [ val(meta), [ bam ] ]
    ch_bai       // channel: [ val(meta), [ bai ] ]
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    ch_fasta_fai // channel: [ val(meta), [ fai ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: freebayes
    //
    ch_bam_bai = ch_bam.join(ch_bai, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, bam, bai -> [meta, bam, bai, [], [], []] }
    FREEBAYES (
        ch_bam_bai,
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
    // MODULE: tabix
    //
    TABIX_TABIX (
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    //
    // MODULE: BCFtools stats
    //
    ch_vcf_tbi = BCFTOOLS_SORT.out.vcf.join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
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

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), [ vcf ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

