#!/usr/bin/env python3

## IDS are organised in directories in the path detailed in os.listdir()
IDS = os.listdir("PAX4_data/")
print(IDS)

rule all:
    input: expand("temp/{id}_Aligned.out.bam", id=IDS)


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
## Run STAR_PE with 6 threads per ID and with 32 threads in total (snakemake -j 32)
rule STAR_PE:
    output:
        # see STAR manual for additional output files
        "temp/{id}_Aligned.out.bam"

    message: "Executing two-pass STAR mapping"
    log:
        "logs/{id}.log"
    threads: 6
    shell:
        """
        # concatenate read 1. Look in the directory containing the samples
        TEST=$(find PAX4_data/{wildcards.id}/*_R1_001.fastq.gz)
        R1=$(echo "${{TEST[*]}}" | paste -sd "," -)

        # concatenate read 2. Look in the directory containing the samples
        TEST=$(find PAX4_data/{wildcards.id}/*_R2_001.fastq.gz)
        R2=$(echo "${{TEST[*]}}" | paste -sd "," -)

        # join both
        PAIRS="$R1 $R2"
        echo $PAIRS

        ## STAR run
        #path to STAR
        /well/mccarthy/production/rna-seq/dependencies/bin/STAR \
        --runThreadN {threads} \
        --genomeDir /well/mccarthy/production/rna-seq/resources/RSEM_GRCh37 \
        --genomeLoad NoSharedMemory \
        --readFilesIn $PAIRS \
        --readFilesCommand zcat \
        --outFileNamePrefix temp/{wildcards.id}_ \
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

rule picard_deduplicate:
    input: expand("temp/{id}_Aligned.out.bam", id=IDS)
    output:
        "bams/{id}.dedup.bam"

    message: "Removing duplicates with PicardTools"
    log:
        "logs/{id}.log"
    threads: 6
    shell:
        """
        # Sort
        samtools sort {input} temp/{wildcards.id}_Aligned.sorted.bam

        ## remove duplicates
        dependencies/picard-tools-2.1.1/picard.jar MarkDuplicates I=temp/{wildcards.id}_Aligned.sorted.bam O={output}

        """
