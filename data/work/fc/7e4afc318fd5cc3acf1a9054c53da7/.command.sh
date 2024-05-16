#!/bin/bash -ue
cutadapt \
    -Z \
    --cores 1 \
     \
    -o SRR18293079_FA-S-P1_1.trim.fastq.gz -p SRR18293079_FA-S-P1_2.trim.fastq.gz \
    SRR18293079_FA-S-P1_1_dummy.fastq SRR18293079_FA-S-P1_2_dummy.fastq \
    > SRR18293079_FA-S-P1.cutadapt.log
cat <<-END_VERSIONS > versions.yml
"CUTADAPT":
    cutadapt: $(cutadapt --version)
END_VERSIONS
