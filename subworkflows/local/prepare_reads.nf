#!/usr/bin/env nextflow
// this subworkflow prepares the inputs from fastq files for the translation



// modules to include in this subworkflow

include {CUTADAPT     } from '/Users/bagordo/Desktop/all/all_bioinformatics/nf-core-nanorepertoire/modules/nf-core/cutadapt/main.nf'
include {FLASH        } from '/Users/bagordo/Desktop/all/all_bioinformatics/nf-core-nanorepertoire/modules/nf-core/flash/main.nf'
include {RENAME       } from '/Users/bagordo/Desktop/all/all_bioinformatics/nf-core-nanorepertoire/modules/local/rename/main.nf'


///Users/bagordo/Desktop/all/all_bioinformatics/nf-core-nanorepertoire/data/*_{1,2}_dummy2.fastq
// main workflow
workflow PREPARE_READS {

    // with take we define the input channels
    take:
    reads      // channel: [mandatory] [ val(meta), [ reads ] ]

    main:
    // define the channels that will be used to store the versions of the software
    input_ch = Channel
            .fromFilePairs(reads)
            .map{ name, reads ->
                [[id:name], [reads[0], reads[1]]]
                }
    ch_versions = Channel.empty()

    CUTADAPT(input_ch)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    FLASH(CUTADAPT.out.reads)
    ch_versions = ch_versions.mix(FLASH.out.versions.first())

    RENAME(FLASH.out.merged, single_end = true)
    ch_versions = ch_versions.mix(RENAME.out.versions.first())


    // defining the output channels that the workflow will
    emit:
    // emitted channels
    trimmed_fastq = CUTADAPT.out.fastq           // channel: [ val(meta), [ fastq ] ]
    merged        = FLASH.out.fastq          // channel: [ val(meta), [ fastq ] ]
    renamed       = RENAMING.out.fastq          // channel: [ val(meta), [ fastq ] ]

    versions      = ch_versions                     // channel: [ versions.yml ]
}

