#!/bin/perl

# original by Martijn van de Bunt
# adapted by Marta Perez Alcantara

##script to compile counts
##should be run on project servers
use strict;
use Getopt::Long;
use Time::Piece;



##default parameters
my $project_dir = (); #data directory variable with default value (NULL)
my $genome = "GRCh37"; #reference genome selection with default value (GRCh37) - can be "GRCh37"|"GRCh38"
my $prefix = (); # prefix for output files with default value (NULL)

#GRCh37 is with GENCODE v19
#GRCh38 is with GENCODE v25

GetOptions("data=s" => \$project_dir, "genome:s" => \$genome,  "prefix:s" => \$prefix);

if(!($project_dir)){ die "No project directory specified - expect name of project directory in data"; }
if($genome ne "GRCh37"){ die "wrong reference genome version specified - can be \"GRCh37\" or \"GRCh38\""; }
my $gencode_file = "../../resources/gencode.v19.annotation.gtf"; #gene annotation file with GENCODE v19 set to default
if(!(-e $project_dir)){ die "Project directory not found - check $project_dir exists"; }

##set flags for processing counts
my $counts = 0;
my $summary = 0;

print "Prefix: $prefix\n";

#get input files
my @count_files = glob("$project_dir/featureCounts/*.gene.counts"); #featureCounts gene counts
my @sum_files = glob("$project_dir/featureCounts/*.gene.counts.summary"); #featureCounts summary files
#set counts to process if files are found
if(scalar(@count_files) > 0){ $counts = 1; }
#and for the summary files
if(scalar(@sum_files) > 0){ $summary = 1; }

#process gene annotation file for gene names
#assumes GENCODE-like annotation
print "Parsing gene annotation file\n";
my %gene_info;
open(my $gencode, "<", $gencode_file);
while(<$gencode>){
    chomp;
    #skip the header
    if($_ =~ m/^#/){ next; }
    my @l = split /\t/;
    #only process the gene entries
    if($l[2] ne "transcript"){ next; }
    #get gene name and biotype
    #currently only using gene name but might be useful to have
    my ($geneid,$genename,$genetype, $transid, $transname) = ();
    if($_ =~ /.*gene_id\s\"([^\"]+).*transcript_id\s\"([^\"]+).*gene_type\s\"([^\"]+).*gene_name\s\"([^\"]+).*transcript_name\s\"([^\"]+).*/){
        $geneid = $1;
        $genename = $4;
        $genetype = $3;
        $transid = $2;
        $transname = $5;
        $transid =~ s/\.\d+$//;
    }
    #ensure the file was parsed correctly
    if(!($geneid) || !($genename) || !($genetype) || !($transid) || !($transname)){ die "GENCODE file parsing error"; }
    #add information to hash
    $gene_info{$geneid}{'name'} = $genename;
    $gene_info{$geneid}{'type'} = $genetype;
    $gene_info{$transid}{'gene'} = $geneid;
    $gene_info{$transid}{'name'} = $transname;
}
close $gencode;

#process count files
if($counts == 1){
    print "Processing gene count files\n";
    #count and size hash
    my (%gene_counts,%sizes,%qc) = ();
    #iterate through the files
    #get gene size information for each gene
    #will use this for TPM calculations
    my @file_header = qw(GeneID GeneName);
    foreach my $count_file (@count_files){
            #get sample name
            my $sample = $count_file;
            $sample =~ s/.*\/(.*).gene.counts/$1/;
            #initialise temp hash for calculating number of reads mapping to the top100 genes
            #useful for QC plot
            my %temp = ();
            push(@file_header,$sample);
            open(my $count, "<", $count_file);
            while(<$count>){
                chomp;
                my @l = split /\t/;
                #skip headers
                if(!($_ =~ m/^ENS/)){ next; }
                #get sizes and ensure they do not differ between files in the same project
                if(exists($sizes{$l[0]}) && $sizes{$l[0]} != $l[5]){ die "Gene size is not equal across gene count files"; }
                $sizes{$l[0]} = $l[5];
                #add the gene count to the gene_counts hash
                $gene_counts{'counts'}{$l[0]}{$sample} = $l[6];
                #add gene count to the temp hash
                $temp{$l[0]} = $l[6];
            }
            close $count;

            #get the first 100 genes and sum across genes
            my $top100 = 0;
            my $sum_all = 0;
            my $n = 1;
            foreach my $gene (sort { $temp{$b} <=> $temp{$a} } keys %temp){
                #get gene counts for the 100 most counted genes
                if($n <= 100){ $top100 += $temp{$gene}; }
                #also get the sum of reads across genes
                $sum_all += $temp{$gene};
                $n++;
            }
            #add the information to the qc hash
            my $top_frac = $top100/$sum_all;
            $top_frac = sprintf("%.3f",$top_frac);
            $qc{$sample}{'top100'} = $top_frac;
    }
    #output the raw counts file
    #also perform first part of TPM calculations
    open(my $count_out, ">", "$prefix.gene.counts.tsv");
    #print new header
    print $count_out join("\t",@file_header)."\n";
    #iterate through the genes
    foreach my $gene (keys %{$gene_counts{'counts'}}){
        print $count_out "$gene\t$gene_info{$gene}{'name'}";
        #run through all samples in the header order
        for(my $i = 2; $i <= $#file_header; $i++){
            #check whether the sample has an entry for that gene
            #then print
            if(exists $gene_counts{'counts'}{$gene}{$file_header[$i]}){
                print $count_out "\t$gene_counts{'counts'}{$gene}{$file_header[$i]}";
                #calculate the RPKs per gene and add up totals per sample
                #RPK = read count per gene divided by gene length in kilobase
                my $rpk = 0;
                #only do so if value isn't zero
                if($gene_counts{'counts'}{$gene}{$file_header[$i]} > 0){
                    $rpk = $gene_counts{'counts'}{$gene}{$file_header[$i]}/($sizes{$gene}/1000);
                }
                #add to hash and keep track of totals per sample
                $gene_counts{'rpks'}{$gene}{$file_header[$i]} = $rpk;
                $gene_counts{'rpk_totals'}{$file_header[$i]} += $rpk;
            }
            #if no entry is found print 0
            else{
                print $count_out "\t0";
                #also create a zero entry for the RPK
                $gene_counts{'rpks'}{$gene}{$file_header[$i]} = 0;
            }
        }
        #print endline
        print $count_out "\n";
    }
    close $count_out;

    #create TPM output file
    open(my $tpm_out, ">", "$prefix.gene.tpm.tsv");
    #print new header
    print $tpm_out join("\t",@file_header)."\n";
    #iterate through the genes
    foreach my $gene (keys %{$gene_counts{'rpks'}}){
        print $tpm_out "$gene\t$gene_info{$gene}{'name'}";
        #run through all samples in the header order
        #have added zeros for all potential missing samples in the previous loop
        for(my $i = 2; $i <= $#file_header; $i++){
            #calculate TPMs, which is RPK divided by the sum of all RPKs per sample divided by 1 million
            #then print
            my $tpm = 0;
            #get the scale factor
            my $scale_factor = $gene_counts{'rpk_totals'}{$file_header[$i]}/1e6;
            #check whether RPK is greater than zero
            if($gene_counts{'rpks'}{$gene}{$file_header[$i]} > 0){
                $tpm = $gene_counts{'rpks'}{$gene}{$file_header[$i]}/$scale_factor;
            }
            #round the TPM to four digits
            $tpm = sprintf("%.4f",$tpm);
            #calculate the number of genes with TPM > 1
            if($tpm > 1){ $qc{$file_header[$i]}{'genes'}++; }
            print $tpm_out "\t$tpm";
        }
        print $tpm_out "\n";
    }
    close $tpm_out;

    #process the summary files for the QC plots
    if($summary != 1){ print "Skipping QC plot generation - no summary files found\n"; }
    else{
        #print the QC hash
        open(my $qc_out, ">", "$prefix.qc.info.tmp");
        print $qc_out "Sample\tfracTop100\texprGenes\n";
        foreach my $smpl (keys %qc){
            print $qc_out "$smpl\t$qc{$smpl}{'top100'}\t$qc{$smpl}{'genes'}\n";
        }
        close $qc_out;
        #hash all the information before printing to a temp file
        my %sum_tmp = ();
        #get the temp header
        my @sum_header = qw(Annotation);
        foreach my $sum_file (@sum_files){
            #get the sample name
            my $name = $sum_file;
            $name =~ s/.*\/(.*).gene.counts.summary/$1/;
            #add it to the temp header
            push(@sum_header,$name);
            #get the information
            open(my $sum_in, "<", $sum_file);
            <$sum_in>;
            while(<$sum_in>){
                chomp;
                my @l = split /\t/;
                #add to the hash
                $sum_tmp{$l[0]}{$name} = $l[1];
            }
            close $sum_in;
        }
        #print the temp output file
        open(my $sum_out, ">", "$prefix.gene.count.qc.summary.tsv");
        print $sum_out join("\t",@sum_header)."\n";
        #print the category information for all samples
        foreach my $category (keys %sum_tmp){
            print $sum_out "$category";
            #print the samples in header order
            for(my $i = 1; $i <= $#sum_header; $i++){
                #check whether sample exists
                #should do for all, but just on the safe side
                if(exists $sum_tmp{$category}{$sum_header[$i]}){ print $sum_out "\t$sum_tmp{$category}{$sum_header[$i]}"; }
                else{ print $sum_out "\t0"; }
            }
            #endline
            print $sum_out "\n";
        }
        close $sum_out;

        #generate the QC plots
        open(my $script, ">", "$prefix.qcplot.R");
        print $script qq^library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
stats <- read.table("$prefix.gene.count.qc.summary.tsv",sep="\\t",header=T,check.names=F)
#remove unused categories
cats <- c("Assigned","Unassigned_MultiMapping","Unassigned_Unmapped","Unassigned_Ambiguity","Unassigned_NoFeatures")
stats <- stats[which(stats[,1] %in% cats),]
#relevel the categories
stats[,1] <- factor(stats[,1],levels=c("Assigned","Unassigned_MultiMapping","Unassigned_Ambiguity","Unassigned_NoFeatures","Unassigned_Unmapped"))
levels(stats[,1]) <- c("Exonic","Multimapping","Ambiguous","NonExonic","Unmapped")
stats_melt <- melt(stats)
p <- ggplot(stats_melt[order(stats_melt[,1]),],aes(x=variable,y=value/1e6,fill=Annotation)) +
geom_bar(stat="identity") +
scale_fill_manual("",breaks=levels(stats[,1]),values=c("#ba4b85","#5d9046","#815ec5","#aa6f35","#cf4a3a")) +
theme_bw() +
ylab("Read pairs (Million)") +
xlab("") +
theme( axis.text.x=element_text(angle=90, vjust=0.5, hjust=1) )
ggsave("$prefix.qc.mapping.category.reads.png",p, width = 17.5, height = 4, units = "in", dpi = 400)


#plot the final QC information
qc <- read.table("$prefix.qc.info.tmp",sep="\\t",header=T,check.names=F)
p <- ggplot(qc,aes(x=fracTop100*100,y=exprGenes)) +
geom_point(cex=3, colour="black") +
geom_text_repel(aes(label=Sample)) +
theme_bw() +
xlab("Reads pairs mapping to top100 genes [%]") +
ylab("Genes with TPM > 1 [#]")
ggsave("$prefix.qc.exprVStop100.png",p)^;
        close $script;
        #run the script
        `Rscript --vanilla $prefix.qcplot.R`;

        #clean up the temp files
        `rm $prefix.qcplot.R $prefix.qc.info.tmp`;
    }
}
else{ print "Skipping count file generation - no input files found in $project_dir\/merged/counts\n"; }

print "Done\n";
