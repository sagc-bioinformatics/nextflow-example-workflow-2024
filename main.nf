params.outdir = 'output'

/*
process MINIMAP2 {

}

process SAMTOOLS_VIEW {

}

process SAMTOOLS_FLAGSTAT {

}

process MULTIQC {

}
*/

process FASTP {

    // we could put this in config also
    // by default, it makes symlinks into the workdir, but this is configurable
    publishDir { "${params.outdir}/QC" }

    // a tuple is just a collection of one or more other values
    // using tuples lets us group metadata and data in a single channel
    input:
    tuple val(id), path(reads)

    // anything we want published, or connected to something downstream,
    // we need to capture in our output
    // it's useful to keep the id so that we can make decisions based on that
    output:
    tuple val(id), path('*.json'), emit: json
    tuple val(id), path('*.html'), emit: html

    // We can pretty much just paste the same block in here
    // nextflow has the same variable substitution syntax so it just works
    """
    # run pre-mapping QC
    echo "Processing pre-mapping QC"
    fastp -i ${reads} \
        --json ${id}.fastp.json \
        --html ${id}.fastp.html
    """
}



workflow {

    // fromFilePairs works too, maybe demonstrate
    // simpleName is because it is a Path, and Groovy lets you convert get*() calls like they're attributes
    // map is an operator - you can use () to call it, but with one arg (a closure), we don't need to
    // Note that you don't need to pass an id around - it just makes it easier!
    channel.fromPath ( 'data/samples/*.fastq' )
        | map { [ it.simpleName, it ] } 
        | FASTP
    
    // FASTP | MINIMAP2 | SAMTOOLS_VIEW | SAMTOOLS_FLAGSTAT | MULTIQC

}