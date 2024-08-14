params.outdir = 'output'
params.genome = 'data/genome.fa'
/*

process MULTIQC {

}
*/

process FASTP {

    // This is a process directive. We could put this in config also
    // by default, it makes symlinks into the workdir, but this is configurable
    // Note that this dir is not controlled in any way! Stuff just gets copied here.
    // I'm also setting a pattern for only publishing html. By default all outputs are copied.
    publishDir "${params.outdir}/QC", pattern: "*.html"

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
    // we can also emit versions here, especially using the new topic channel syntax might be nice to show...

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

process MINIMAP2 {

    // we won't have publishDir here, because we don't care about sam files
    cpus 2 // we can set this in config also, but it will mean the implicit variable task.cpus is set

    input:
    tuple val(id), path(reads)
    path genome // a separate channel, not associated with a particular sample

    output:
    tuple val(id), path('*.sam'), emit: sam

    // again we can just paste our script verbatim
    // Nextflow keeps info on what commands are run anyway, so we can skip that bit (show how to inspect)
    """
    # map reads
    echo "mapping reads"
    minimap2 -t ${task.cpus} -a -x sr ${genome} ${reads} > ${id}.sam 2> ${id}.minimap2.log
    """
}

process SAMTOOLS_VIEW {
    
    publishDir { "${params.outdir}/aligned" }

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path('*.bam'), emit: bam

    """
    # convert SAM to BAM
    echo "convert SAM to BAM"
    samtools view -bh -o ${id}.bam ${sam}
    """
}

process SAMTOOLS_FLAGSTAT {

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('*.flagstat.txt')

    """
    # post-mapping QC
    samtools flagstat ${bam} > ${id}.flagstat.txt
    """
}

workflow {

    //
    // SET UP CHANNELS
    //

    // fromFilePairs works too, maybe demonstrate
    // simpleName is because it is a Path, and Groovy lets you convert get*() calls like they're attributes
    // map is an operator - you can use () to call it, but with one arg (a closure), we don't need to
    // Note that you don't need to pass an id around - it just makes it easier!
    ch_input = channel.fromPath ( 'data/samples/*.fastq' )
        | map { [ it.simpleName, it ] }

    // you could use params.genome in a process without making a channel
    // but there are problems with that approach
    // I've done channel.value explicitly here, but you can omit that - it happens implicitly if you pass a value
    ch_genome = channel.value( file(params.genome, checkIfExists: true) )

    //
    // CONNECT PROCESSES
    //

    FASTP ( ch_input )
    
    MINIMAP2 ( ch_input, ch_genome ) | SAMTOOLS_VIEW | SAMTOOLS_FLAGSTAT

    // TODO: finally, multiQC
    // Should I demonstrate the use of topic channels?

}