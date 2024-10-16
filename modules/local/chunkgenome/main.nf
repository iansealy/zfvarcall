process CHUNK_GENOME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0':
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(fai)  // channel: [ val(meta), [ fai ] ]

    output:
    path "*.bed"        , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    awk ${args} -vFS="\t" '{
        contig_size = \$2
        if (name == "" || (chunk_size > ${params.split_genome})) {
            # Start a new chunk
            name = sprintf("%d.bed", ++i)
            chunk_size = 0
        }
        chunk_size += contig_size
        print \$1 "\t0\t" \$2 > name
    }' ${fai}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch 1.bed 2.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
