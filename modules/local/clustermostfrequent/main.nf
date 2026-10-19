process CLUSTER_MOSTFREQUENT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.78'
        : 'biocontainers/biopython:1.78'}"

    input:
    tuple val(meta), path(sequences), path(clusters)

    output:
    tuple val(meta), path("*_most_frequent.fasta"), emit: fasta
    tuple val(meta), path("*_most_frequent.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cluster_most_frequent.py \\
        ${sequences} \\
        ${clusters} \\
        ${prefix}_most_frequent \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
