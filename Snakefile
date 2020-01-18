#!/usr/bin/env python3
import os




## IDS are organised in directories in the path detailed in os.listdir()
IDS = os.listdir("PAX4_data/")
print(IDS)

rule all:
    input: "jan_2020.gene.counts.tsv"


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


## add to bashrc the following if picardtools runs out of memory
#export _JAVA_OPTIONS="-Xms4g -Xmx4g -XX:ParallelGCThreads=1"

rule picard_deduplicate:
    input: "temp/{id}_Aligned.out.bam"
    output:
        dedup="bams/{id}.dedup.bam",
        sort="temp/{id}_Aligned.sorted.bam",
        metrics="temp/{id}_marked_dup_metrics.txt"

    message: "Removing duplicates with PicardTools. Make sure to give enough memory to java with export _JAVA_OPTIONS= [comma] -Xms4g -Xmx4g -XX:ParallelGCThreads=1 [comma]. Run on cluster with snakemake  --jobs 164 --cluster [comma] qsub -cwd -V -pe shmem 10 -P mccarthy.prjc -q short.qc -N picard [comma]"
    log:
        "logs/{id}_dedup.log"
    shell:
        """
        # Sort
        samtools sort {input} -o temp/{wildcards.id}_Aligned.sorted.bam &> {log}

        ## remove duplicates
        java -jar /well/mccarthy/users/martac/bin/picard-tools-2.1.1/picard.jar MarkDuplicates I={output.sort} O={output.dedup} M={output.metrics} &>> {log}

        ## cleanup

        #rm picard.e*
        #rm picard.o*
        """

rule counts_featureCounts:
    input: "bams/{id}.dedup.bam"
    output:
        "featureCounts/{id}.gene.counts"

    message: "Calculate raw counts per gene with featureCounts. Run on cluster with snakemake  --jobs 164 --cluster [comma] qsub -cwd -V -pe shmem 6 -P mccarthy.prjc -q short.qc -N featurecounts [comma]"
    log:
        "logs/{id}_featureCounts.log"
    threads: 6
    shell:
        """
        ## indexing for featureCounts
        samtools index {input}

        /well/mccarthy/production/rna-seq/dependencies/bin/featureCounts -p -t exon -g gene_id -s 0 -T 4 -B -C -a /well/mccarthy/production/rna-seq/resources/gencode.v19.annotation.gtf -o {output} {input} &> {log}

        # featurecounts.e*
        #  featurecounts.o*
        """

rule tidy_counts_TPM:
    input: expand("featureCounts/{id}.gene.counts",id=IDS)
    output:
        "jan_2020.gene.counts.tsv"

    message: "Create transcript per million (TPM) count table. Also produces tidy raw count table for downstream analyses and QC plots."
    log:
        "logs/tidy.log"
    threads: 6
    shell:
        """
        module load R/3.2.2
        ## running perl script
        /apps/well/perl/5.16.3/bin/perl scripts/TPM_qc.pl --data ./ --genome GRCh37 --prefix jan_2020 &> {log}

        """
