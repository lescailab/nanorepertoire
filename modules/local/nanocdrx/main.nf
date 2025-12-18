process NANOCDRX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/nanocdr-x:1.0.0'
        : 'docker.io/lescailab/nanocdr-x:1.0.0'}"

    input:
    tuple val(meta), path(translated)

    output:
    tuple val("${meta.id}"), val("${meta.immunisation}"), val("${meta.boost}"), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.hist"), emit: hist
    tuple val(meta), path("*.tsv"), emit: tsv
    path '*.hist', emit: histonly
    path '*.tsv', emit: tsvonly
    val meta, emit: metaonly
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_nanocdrx.py \\
        -i ${translated} \\
        -p ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocdr-x: \$(predict_cdrs --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """
}
