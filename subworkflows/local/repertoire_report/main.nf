//
// Subworkflow: REPERTOIRE_REPORT
//
include { AGGREGATE_STATS      } from '../../../modules/local/aggregate_stats/main'
include { NANOREPERTOIRE_REPORT } from '../../../modules/local/nanorepertoire_report/main'

workflow REPERTOIRE_REPORT {

    take:
    report_template     // path: analysis_report.qmd
    loop_tree_template  // path: loop_tree.qmd
    clusteread          // channel: [ val(meta), path(summary) ]
    cdrhistograms       // channel: [ val(meta), path(hist) ]
    cdrtsv              // channel: [ val(meta), path(tsv) ]
    metadata            // channel: [ val(meta) ]

    main:

    ch_versions = Channel.empty()

    AGGREGATE_STATS (
        report_template,
        loop_tree_template,
        clusteread,
        cdrhistograms,
        cdrtsv,
        metadata
    )
    ch_versions = ch_versions.mix(AGGREGATE_STATS.out.versions)

    NANOREPERTOIRE_REPORT (
        AGGREGATE_STATS.out.clustercounts,
        AGGREGATE_STATS.out.cdrcounts,
        AGGREGATE_STATS.out.cdrhists,
        AGGREGATE_STATS.out.clusterbig,
        AGGREGATE_STATS.out.fastaseq
    )
    ch_versions = ch_versions.mix(NANOREPERTOIRE_REPORT.out.versions)

    emit:
    report   = NANOREPERTOIRE_REPORT.out.report  // path: *.html
    versions = ch_versions                       // channel: [ versions.yml ]
}
