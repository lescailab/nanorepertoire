process REPORTSAMPLES {
    label 'process_single'

    input:
    val(metadata)

    output:
    path("sampledata.tsv"), emit: reportsamples

    exec:
    def sampleData = new File("${task.workDir}/sampledata.tsv")
    sampleData.append('ID\tindividual\timmunisation\tboost\n')

    // Popola il file sampledata.tsv con i dati dalla mappa metadata
    metadata.each() { map ->
        sampleData.append("${map.sampleID}\t${map.individualID}\t${map.immunisation}\t${map.boost}\n")
    }
}
