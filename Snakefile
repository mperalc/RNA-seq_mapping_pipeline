#!/usr/bin/env python3

## join ls output by comma using bash
#TEST=$(ls test/RHP7393/*.gz)
#bar=$(printf "{,}%s" "${TEST[@]}")
#echo ${TEST[*]// /|}
#echo $bar

#TEST=$(find test/RHP7393/*.gz)
#bar=$(printf "{,}%s" "${TEST[@]}")
#echo $bar
IDS = "RHP7393 RHP7394".split()
print(IDS)

rule all:
    input:
        expand("test/{id}/{id}_output.tsv", id=IDS)

rule test_wildcards:
    output: "test/{id}/{id}_output.tsv"
    shell:
        """
        printf '%s\n' test/{wildcards.id}/*_R1_001.fastq.gz  ##verify path, it's very important
        TEST=$(find test/{wildcards.id}/*_R1_001.fastq.gz)
        bar=$(printf "{{,}}%s" "${{TEST[@]}}")
        echo $bar
        touch {output}
        """

rule STAR_PE:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = ["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = ["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"] #optional
    output:
        # see STAR manual for additional output files
        "star/pe/{sample}/Aligned.out.sam"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index="index",
        # optional parameters
        extra=""
    threads: 8
    shell:
        """
        #path to STAR
        dependencies/bin/STAR  \
        #number of threads
        --runThreadN 6 \
        #path to directory with genome indices
        --genomeDir resources/RSEM_GRCh37 \
        #shared memory is not used to load the genome
        --genomeLoad NoSharedMemory \
        #path to files containing sequences to be mapped
        --readFilesIn data/test/fastq/SRR1027171_1.fastq.gz    data/test/fastq/SRR1027171_2.fastq.gz \
        #command to decompress fastq files
        --readFilesCommand zcat \
        #prefix of output file names
        --outFileNamePrefix  data/test/temp/SRR1027171/SRR1027171 \
        #outputs multiple alignments for each read in random order and randomizes the primary alignment from the highest scoring alignments
        --outMultimapperOrder Random  \
        #output unsorted Aligned.out.bam ﬁle
        --outSAMtype BAM   Unsorted \
        # the standard SAM attributes
        --outSAMattributes NH   HI   AS   NM   MD \
        #output unmapped reads within the main SAM ﬁle (i.e. Aligned.out.sam)
        --outSAMunmapped Within \
        # the MAPQ value for unique mappers (default = 255)
        --outSAMmapqUnique 50 \
        #  SAM/BAM read group line
        --outSAMattrRGline ID:SRR1027171   SM:sample   CN:test   PL:ILLUMINA \
        --outSAMheaderHD @HD   VN:1.4   SO:unsorted  \
        # max number of multiple alignments for a read that will be output to the SAM/BAM ﬁles.
        --outSAMmultNmax 1 \
        #ENCODE standard: reduces the number of spurious junctions
        --outFilterType BySJout \
        #ENCODE standard: max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
        --outFilterMultimapNmax 20 \
        #ENCODE standard: maximum number of mismatches per pair, large number switches oﬀ this ﬁlter
        --outFilterMismatchNmax 999 \
        #ENCODE standard: max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.04*200=8 for the paired read
        --outFilterMismatchNoverLmax 0.04 \
        #ENCODE standard: min intron length
        --alignIntronMin 20 \
        #ENCODE standard: max intron length
        --alignIntronMax 1000000 \
        #ENCODE standard: max genomic distance between mates
        --alignMatesGapMax 1000000 \
        #ENCODE standard: minimum overhang for unannotated junctions
        --alignSJoverhangMin 8 \
        #ENCODE standard: minimum overhang for annotated junctions
        --alignSJDBoverhangMin 1 \
        #extra alignment score for alignmets that cross database junctions
        --sjdbScore 1 \
        #basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the ﬂy
        --twopassMode Basic
        """
