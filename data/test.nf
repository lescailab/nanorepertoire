// Cutadapt finds and removes adapter sequences, primers, 
// poly-A tails and other types of unwanted sequence from 
// your high-throughput sequencing reads. 



input_ch = Channel
.fromFilePairs('/Users/bagordo/Desktop/all_bioinformatics/nf-core-nanorepertoire/data/*_{1,2}_dummy.fastq')
.map{ name, reads ->
    [[id:name], [reads[0], reads[1]]]
}

include {CUTADAPT} from '/Users/bagordo/Desktop/all_bioinformatics/nf-core-nanorepertoire/modules/nf-core/cutadapt/main.nf'

input_ch.dump(tag: 'input_ch')

workflow {
    letters_ch = CUTADAPT(input_ch)
}