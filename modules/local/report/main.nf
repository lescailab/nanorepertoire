process REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-genomicfeatures_bioconductor-ggtree_bioconductor-gviz_pruned:d368c94e0967e319':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-genomicfeatures_bioconductor-ggtree_bioconductor-gviz_pruned:fdfe24a734a90d8f' }"

    input:
    path (report)
    path (loop_tree)
    path (clustersummaries)
    path (cdr3histograms)
    path (cdr3tables)
    val  (metadata)
    val  (options)

    output:
    path "*.RData", emit: rdata
    path "*.html", emit: report
    path "cdr3_boost_overview_table.tsv", emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sampleData = new File("${workDir}/sampledata.tsv")
    sampleData.append("ID\tindividual\timmunisation\tboost\n")

    // Popola il file sampledata.tsv con i dati dalla mappa metadata
    metadata.each() { map ->
        sampleData.append("${map.sampleID}\t${map.individualID}\t${map.immunisation}\t${map.boost}\n")
    }

    // Concatena le liste dei file in stringhe separate da virgola
    clusterList = clustersummaries.join(",")
    histoList = cdr3histograms.join(",")
    tableList = cdr3tables.join(",")

    def output = report.baseName

    """
    ln -s '$moduleDir/nibsc_report.css' .
    ln -s '${workDir}/sampledata.tsv' .

    quarto render analysis_report.qmd \\
        -P clusterList:\"$clusterList\" \\
        -P histoList:\"$histoList\" \\
        -P tableList:\"$tableList\" \\
        -P sampleData:\"sampledata.tsv\" \\
        -P sizeThreshold:\"${params.cluster_size_threshold}\" \\
        -P loopFile:\"loop_tree.qmd" \\
        -P calcTree:\"${params.calculate_tree}\" \\
        --output "${output}.html"
    """

}
