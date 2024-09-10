// This subworkflows is used to create the report of nanorepertoire analysis
// modules to include in this subworkflow

include { MULTIQC as MULTIQC1   } from '../../../modules/nf-core/multiqc/main.nf'
include { MULTIQC as MULTIQC2   } from '../../../modules/nf-core/multiqc/main.nf'
include { REPORT    }             from '../../../modules/local/report/main.nf'

workflow NANOREPORT {

    // with take we define the input channels

    take:
    fastqc_zip
    cutadapt_log
    report
    loop_tree
    clustersummaries
    cdr3histograms
    cdr3tables
    metadata
    options


    main:

    ch_versions = Channel.empty()

    MULTIQC1(fastqc_zip,   [], [],[], [], [])
    ch_versions = ch_versions.mix(MULTIQC1.out.versions.first())

    MULTIQC2(cutadapt_log, [], [],[], [], [])
    ch_versions = ch_versions.mix(MULTIQC2.out.versions.first())

    REPORT(report, loop_tree, clustersummaries, cdr3histograms, cdr3tables, metadata, options)



    emit:
    // TODO nf-core: edit emitted channels
    multiqc_report     = MULTIQC1.out.report            // channel: [ val(meta), [ bam ] ]
    multiqc_report1    = MULTIQC2.out.report            // channel: [ val(meta), [ bam ] ]
    analysis_report    = REPORT.out.report              // channel: [ val(meta), [ bai ] ]

    versions           = ch_versions                    // channel: [ versions.yml ]
}

