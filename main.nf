process MINIMAP2 {

}

process SAMTOOLS_VIEW {

}

process SAMTOOLS_FLAGSTAT {

}

process FASTP {

}

process MULTIQC {

}

workflow {

    FASTP | MINIMAP2 | SAMTOOLS_VIEW | SAMTOOLS_FLAGSTAT | MULTIQC

}