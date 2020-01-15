#!/usr/bin/env python3
import os

## IDS are organised in directories in the path detailed in os.listdir()
IDS = os.listdir("PAX4_data/")
print(IDS)

rule all:
    input: expand("bams/{id}.dedup.bam", id=IDS)


## Run STAR_PE with 6 threads per ID and with 36 threads in total (snakemake -j 36)
rule STAR_PE:
    output:
        # see STAR manual for additional output files
        "temp/{id}_Aligned.out.bam"

    message: "Executing two-pass STAR mapping. Run on cluster with snakemake  --jobs 32 --cluster [comma] qsub {params.environment} {params.project} {params.queue} {params.name}[comma]"
    log:
        "logs/{id}_map.log"
    threads: 6
    params:
        environment= "-cwd -V -pe shmem 6",
        project= "-P mccarthy.prjc",
        queue="-q long.qc",
        name= "-N STAR_PE"
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
        --twopassMode Basic &> {log}

        ## cleanup
        rm STAR_PE.e*
        rm STAR_PE.o*

        """



rule picard_deduplicate:
    input: "temp/{id}_Aligned.out.bam"
    output:
        "bams/{id}.dedup.bam"

    message: "Removing duplicates with PicardTools. Run on cluster with snakemake  --jobs 1 --cluster [comma] qsub -cwd -V -pe shmem 6 -P mccarthy.prjc -q long.qc -N picard [comma]"
    log:
        "logs/{id}_dedup.log"
    threads: 6
    shell:
        """
        # Sort
        samtools sort {input} temp/{wildcards.id}_Aligned.sorted.bam

        ## remove duplicates
        dependencies/picard-tools-2.1.1/picard.jar MarkDuplicates I=temp/{wildcards.id}_Aligned.sorted.bam O={output} &> {log}

        """

rule counts_featureCounts:
    input: "bams/{id}.dedup.bam"
    output:
        "featureCounts/{id}.gene.counts"

    message: "Calculate raw counts per gene with featureCounts."
    log:
        "logs/{id}_featureCounts.log"
    threads: 6
    shell:
        """
        ## indexing for featureCounts
        samtools index {input}

        dependencies/bin/featureCounts -p -t exon -g gene_id -s 0 -T 4 -B -C -a resources/gencode.v19.annotation.gtf -o {output} {input} &> {log}

        """

rule tidy_counts_TPM:
    input: "featureCounts/{id}.gene.counts"
    output:
        "counts/{id}.gene.counts.tsv"

    message: "Create transcript per million (TPM) count table. Also produces tidy raw count table for downstream analyses and QC plots."
    log:
        "logs/{id}_tidy.log"
    threads: 6
    shell:
        """
        ## running perl script
        perl TPM_qc.pl –data  –genome GRCh37

        dependencies/bin/featureCounts -p -t exon -g gene_id -s 0 -T 4 -B -C -a resources/gencode.v19.annotation.gtf -o {output} {input} &> {log}

        """
