// this subworkflow use fasta file input to analyse and study cd loop of nanobodies
// modules to include in this subworkflow

include { CDHIT_CDHIT } from '../../../modules/nf-core/cdhit/cdhit/main.nf'
include { READCDHIT } from '../../../modules/local/readcdhit/main.nf'
include { NANOCDRX } from '../../../modules/local/nanocdrx/main.nf'
//include {MAFFT         } from '../../../modules/nf-core/mafft/main.nf'

// main workflow
workflow FASTA_CLUSTERING {
    take:
    input // channel: [[id], [reads_forward, reads_reverse]]

    main:
    // define the channels that will be used to store the versions of the software

    ch_versions = Channel.empty()
    input_ch = input


    CDHIT_CDHIT(input_ch)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions.first())

    READCDHIT(CDHIT_CDHIT.out.clusters)
    ch_versions = ch_versions.mix(READCDHIT.out.versions.first())

    NANOCDRX(CDHIT_CDHIT.out.fasta)
    ch_versions = ch_versions.mix(NANOCDRX.out.versions.first())

    emit:
    cluster = CDHIT_CDHIT.out.fasta // channel: [[id], [ fastq_merged ]]
    clusteread = READCDHIT.out.summaryonly // channel: [[id], [ fastq_renamed]]
    cdrpredicted = NANOCDRX.out.fasta
    cdrtsv = NANOCDRX.out.tsvonly
    cdrhistograms = NANOCDRX.out.histonly
    cdrmeta = NANOCDRX.out.metaonly
    versions = ch_versions // channel: [ versions.yml ]
}
