//
// Subworkflow to prepare reference genome files
//

include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { BWA_INDEX                      } from '../../modules/nf-core/bwa/index/main'

workflow PREPARE_GENOME {

    take:
    fasta                                       // string: /path/to/genome.fasta
    fasta_fai                                   // string: /path/to/genome.fasta.fai
    fasta_dict                                  // string: /path/to/genome.dict
    bwa_index                                   // string: /path/to/bwa.index.dir

    main:

    fasta = file(fasta, checkIfExists: true)
    ch_fasta = [[id:fasta.baseName], fasta]

    SAMTOOLS_FAIDX (
        ch_fasta,
        [ [ id:'fai' ], [] ]
    )

    GATK4_CREATESEQUENCEDICTIONARY (
        ch_fasta
    )

    BWA_INDEX (
        ch_fasta
    )

    emit:
    fasta      = ch_fasta                                 // channel: [ val(meta), [ fasta ] ]
    fasta_fai  = SAMTOOLS_FAIDX.out.fai                   // channel: [ val(meta), [ fai ] ]
    fasta_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict  // channel: [ val(meta), [ dict ] ]
    bwa_index  = BWA_INDEX.out.index                      // channel: [ val(meta), [ bwa ] ]
}

