/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { SEQKIT_SPLIT2          } from '../modules/nf-core/seqkit/split2/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SPLITS } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_LANES  } from '../modules/nf-core/samtools/merge/main'
include { GATK4_ADDORREPLACEREADGROUPS            } from '../modules/nf-core/gatk4/addorreplacereadgroups/main'
include { BIOBAMBAM_BAMSORMADUP  } from '../modules/nf-core/biobambam/bamsormadup/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { GATK4_HAPLOTYPECALLER  } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_zfvarcall_pipeline'

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
    // MODULE: fastp
    //
    FASTP (
        ch_samplesheet,
        [],    // no need to specify adapters
        false, // no need to discard passing reads
        false, // no need to keep failed reads
        false  // no need to merge reads
    )
    ch_reads_trimmed = FASTP.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC (
        ch_reads_trimmed
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: SeqKit split2
    //
    SEQKIT_SPLIT2 (
        ch_reads_trimmed
    )
    ch_reads_split = SEQKIT_SPLIT2.out.reads.transpose()
        .map{ meta, read -> [meta + [split: get_split_num(read.getName())], read]}
        .groupTuple()
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    //
    // MODULE: BWA mem
    //
    BWA_MEM (
        ch_reads_split,
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

    //
    // MODULE: biobambam2 bamsormadup
    //
    BIOBAMBAM_BAMSORMADUP (
        GATK4_ADDORREPLACEREADGROUPS.out.bam,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BIOBAMBAM_BAMSORMADUP.out.metrics.collect{it[1]})
    ch_versions = ch_versions.mix(BIOBAMBAM_BAMSORMADUP.out.versions)

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

    //
    // MODULE: mosdepth
    //
    ch_bam_bai_1 = SAMTOOLS_MERGE_LANES.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
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
    ch_bam_bai_0 = SAMTOOLS_MERGE_LANES.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, bam, bai -> [meta, bam, bai] }
    SAMTOOLS_FLAGSTAT (
        ch_bam_bai_0
    )
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    // TODO: Use intervals
    //
    // MODULE: GATK HaplotypeCaller
    //
    ch_bam_bai_2 = SAMTOOLS_MERGE_LANES.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, bam, bai -> [meta, bam, bai, [], []] }
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_2,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        [[], []], // no need for dbSNP file
        [[], []]  // no need for dbSNP index
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

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

// Helper functions

// Extract split number from a read filename
// e.g. Get "2" from "ERS01_1.fastp.part_002.fastq.gz"
def get_split_num(String read) {
    return (read =~ /(\d+)\.(fastq|fq)(\.gz)?$/)[0][1].toInteger()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
