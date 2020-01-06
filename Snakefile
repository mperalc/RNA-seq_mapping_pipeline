#!/usr/bin/env python3

IDS = "RHP7393 RHP7394".split()
print(IDS)

rule all:
    input: expand("test/{id}/Aligned.out.sam", id=IDS)


rule test_wildcards:
    output: "test/{id}/{id}_output.tsv"
    shell:
        """
        printf '%s\n' test/{wildcards.id}/*_R1_001.fastq.gz  ##verify path, it's very important
        # concatenate read 1
        TEST=$(find test/{wildcards.id}/*_R1_001.fastq.gz)
        R1=$(echo "${{TEST[*]}}" | paste -sd "," -)
        echo $R1

        # concatenate read 2
        TEST=$(find test/{wildcards.id}/*_R2_001.fastq.gz)
        R2=$(echo "${{TEST[*]}}" | paste -sd "," -)
        echo $R2

        # join both
        PAIRS="$R1 $R2"
        echo $PAIRS
        echo $PAIRS > test/test.txt

        touch {output}
        """

rule STAR_PE:
    # input:
    #     # use a list for multiple fastq files for one sample
    #     # usually technical replicates across lanes/flowcells
    #     fq1 = ["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
    #     # paired end reads needs to be ordered so each item in the two lists match
    #     fq2 = ["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"] #optional
    output:
        # see STAR manual for additional output files
        "test/{id}/Aligned.out.sam"

    message: "Executing two-pass STAR mapping"
    log:
        "test/logs/{id}.log"
    threads: 8
    shell:
        """
        # concatenate read 1
        TEST=$(find test/{wildcards.id}/*_R1_001.fastq.gz)
        R1=$(echo "${{TEST[*]}}" | paste -sd "," -)

        # concatenate read 2
        TEST=$(find test/{wildcards.id}/*_R2_001.fastq.gz)
        R2=$(echo "${{TEST[*]}}" | paste -sd "," -)

        # join both
        PAIRS="$R1 $R2"
        echo $PAIRS

        ## STAR run
        #path to STAR
        /well/mccarthy/production/rna-seq/dependencies/bin/STAR \
        --runThreadN 6 \
        --genomeDir /well/mccarthy/production/rna-seq/resources/RSEM_GRCh37 \
        --genomeLoad NoSharedMemory \
        --readFilesIn $PAIRS \
        --readFilesCommand zcat \
        --outFileNamePrefix test/{wildcards.id}/temp \
        --outMultimapperOrder Random \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI AS NM MD \
        --outSAMunmapped Within \
        --outSAMmapqUnique 50 \
        --outSAMattrRGline ID:{wildcards.id} SM:sample CN:test PL:ILLUMINA \
        --outSAMheaderHD @HD VN:1.4 SO:unsorted \
        --outSAMmultNmax 1 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --twopassMode Basic

        """
