// this subworkflow prepares the inputs from fastq files for the translation
// modules to include in this subworkflow

include {CDHIT_CDHIT   } from '../../../modules/nf-core/cdhit/cdhit/main.nf'
include {READCDHIT     } from '../../../modules/local/readcdhit/main.nf'
include {GETCDR3       } from '../../../modules/local/getcdr3/main.nf'


///Users/bagordo/Desktop/all/all_bioinformatics/nf-core-nanorepertoire/data/*_{1,2}_dummy2.fastq
// main workflow
workflow FASTQTOFASTA_CLUSTERING_CDR3 {

    // with take we define the input channels
    take:
    input      // channel: [[id], [reads_forward, reads_reverse]]
    single_end

    main:
    // define the channels that will be used to store the versions of the software

    ch_versions = Channel.empty()
    input_ch = input


    CDHIT_CDHIT(input_ch)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions.first())

    READCDHIT(CDHIT_CDHIT.out.clusters, single_end = true)
    ch_versions = ch_versions.mix(READCDHIT.out.versions.first())

    GETCDR3(CDHIT_CDHIT.out.fasta, single_end = true)
    ch_versions = ch_versions.mix(GETCDR3.out.versions.first())




    // defining the output channels that the workflow will
    emit:
    // emitted channels
    cluster          = CDHIT_CDHIT.out.fasta              // channel: [[id], [ fastq_merged ]]
    clusteread       = READCDHIT.out.summary           // channel: [[id], [ fastq_renamed]]
    cdrpredicted     = GETCDR3.out.fasta
    cdrpredicted     = GETCDR3.out.tsv
   
    versions      = ch_versions                     // channel: [ versions.yml ]
}

