/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_FASTP_FASTQC_SPLIT        } from '../subworkflows/local/fastq_fastp_fastqc_split'
include { FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX } from '../subworkflows/local/fastq_align_bwa_merge_addrg_markdup_merge_index'
include { MOSDEPTH                        } from '../modules/nf-core/mosdepth/main'
include { SAMTOOLS_FLAGSTAT               } from '../modules/nf-core/samtools/flagstat/main'
include { BAM_GATK4_HAPLOTYPECALLER_STATS } from '../subworkflows/local/bam_gatk4_haplotypecaller_stats'
include { BAM_FREEBAYES_SORT_INDEX_STATS  } from '../subworkflows/local/bam_freebayes_sort_index_stats'
include { BCFTOOLS_MPILEUP                } from '../modules/nf-core/bcftools/mpileup/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_zfvarcall_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ZFVARCALL {

    take:
    ch_samplesheet        // channel: samplesheet read in from --input
    ch_fasta              // channel: [ val(meta), [ fasta ] ]
    ch_fasta_fai          // channel: [ val(meta), [ fai ] ]
    ch_fasta_dict         // channel: [ val(meta), [ dict ] ]
    ch_bwa_index          // channel: [ val(meta), [ bwa ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: fastp, FastQC and SeqKit split2
    //
    FASTQ_FASTP_FASTQC_SPLIT (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(FASTQ_FASTP_FASTQC_SPLIT.out.versions)

    //
    // SUBWORKFLOW: BWA mem, Samtools merge (splits), GATK AddOrReplaceReadGroups,
    // biobambam2 bamsormadup, Samtools merge (lanes) and Samtools index
    //
    FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX (
        FASTQ_FASTP_FASTQC_SPLIT.out.reads,
        ch_bwa_index,
        ch_fasta,
        ch_fasta_fai
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.multiqc.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.versions)

    //
    // MODULE: mosdepth
    //
    ch_bam_bai_1 = FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bam.join(FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, bam, bai -> [meta, bam, bai, []] }
    MOSDEPTH (
        ch_bam_bai_1,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    //
    // MODULE: Samtools flagstat
    //
    ch_bam_bai_0 = FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bam.join(FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, bam, bai -> [meta, bam, bai] }
    SAMTOOLS_FLAGSTAT (
        ch_bam_bai_0
    )
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    //
    // SUBWORKFLOW: GATK HaplotypeCaller and BCFtools stats
    //
    BAM_GATK4_HAPLOTYPECALLER_STATS (
        FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bam,
        FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bai,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict
    )
    ch_versions = ch_versions.mix(BAM_GATK4_HAPLOTYPECALLER_STATS.out.versions)

    //
    // SUBWORKFLOW: freebayes, VCF sort, index and stats
    //
    BAM_FREEBAYES_SORT_INDEX_STATS (
        FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bam,
        FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bai,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(BAM_FREEBAYES_SORT_INDEX_STATS.out.versions)

    //
    // MODULE: BCFtools mpileup
    //
    ch_bam_1 = FASTQ_ALIGN_BWA_MERGE_ADDRG_MARKDUP_MERGE_INDEX.out.bam.map{ meta, bam -> [meta, bam, []] }
    BCFTOOLS_MPILEUP (
        ch_bam_1,
        ch_fasta,
        false // no need to save mpileup
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_replace_names = params.multiqc_replace_names ?
        Channel.fromPath(params.multiqc_replace_names, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_sample_names  = params.multiqc_sample_names ?
        Channel.fromPath(params.multiqc_sample_names, checkIfExists: true) :
        Channel.empty()


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_multiqc_replace_names.toList(),
        ch_multiqc_sample_names.toList()
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
