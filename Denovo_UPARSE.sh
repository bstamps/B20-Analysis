#!/bin/bash

#To analyze multiple Illumina MiSeq runs using the most recent primer set, a slightly different pipeline is required
#This includes analyzing the 16S and 18S amplified together independently. 

# Use QIIME to prepare your FASTQ files. No reverse primer is required for 18S reads. 
extract_barcodes.py -f seq/splitRuns/Run1_Joined.fastq  -a -m tags.txt -l 12 -o seq/splitRuns/preppedRun1/
extract_barcodes.py -f seq/splitRuns/Run2_Joined.fastq  -a -m tags.txt -l 12 -o seq/splitRuns/preppedRun2/
extract_barcodes.py -f seq/splitRuns/Run1_Euk1.fastq  -l 12 -o seq/splitRuns/preppedRun1EukA/
extract_barcodes.py -f seq/splitRuns/Run1_Euk2.fastq  -l 12 -o seq/splitRuns/preppedRun1EukB/
extract_barcodes.py -f seq/splitRuns/Run2_Euk1.fastq  -l 12 -o seq/splitRuns/preppedRun2EukA/
extract_barcodes.py -f seq/splitRuns/Run2_Euk2.fastq  -l 12 -o seq/splitRuns/preppedRun2EukB/

mkdir seq/splitRuns/preppedRun1Euk
mkdir seq/splitRuns/preppedRun2Euk
mkdir seq/SlOut/
mkdir seq/SlOutEuk/
cat seq/splitRuns/preppedRun1EukA/reads.fastq seq/splitRuns/preppedRun1EukB/reads.fastq > seq/splitRuns/preppedRun1Euk/reads.fastq
cat seq/splitRuns/preppedRun1EukA/barcodes.fastq seq/splitRuns/preppedRun1EukB/barcodes.fastq > seq/splitRuns/preppedRun1Euk/barcodes.fastq
cat seq/splitRuns/preppedRun2EukA/reads.fastq seq/splitRuns/preppedRun2EukB/reads.fastq > seq/splitRuns/preppedRun2Euk/reads.fastq
cat seq/splitRuns/preppedRun2EukA/barcodes.fastq seq/splitRuns/preppedRun2EukB/barcodes.fastq > seq/splitRuns/preppedRun2Euk/barcodes.fastq

# Use QIIME to demultiplex the data, with -q 0. Store output as fastq format (we will quality filter with usearch). This part of the analysis was taken from suggestions within the QIIME forums on best practices to integrate UPARSE into a usable QIIME workflow.
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun1/reads.fastq -b seq/splitRuns/preppedRun1/barcodes.fastq -m Run1Tags.txt --barcode_type 12 -o seq/splitRuns/SlOutRun1/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun2/reads.fastq -b seq/splitRuns/preppedRun2/barcodes.fastq -m Run2Tags.txt --barcode_type 12 -o seq/splitRuns/SlOutRun2/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun1Euk/reads.fastq -b seq/splitRuns/preppedRun1Euk/barcodes.fastq -m Run1Tags.txt --barcode_type 12 -o seq/splitRuns/SlOutRun1Euk/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun2Euk/reads.fastq -b seq/splitRuns/preppedRun2Euk/barcodes.fastq -m Run2Tags.txt --barcode_type 12 -o seq/splitRuns/SlOutRun2Euk/

cat seq/splitRuns/SlOutRun1/seqs.fna seq/splitRuns/SlOutRun2/seqs.fna > seq/SlOut/seqs.fna
cat seq/splitRuns/SlOutRun1/seqs.fastq seq/splitRuns/SlOutRun2/seqs.fastq > seq/SlOut/seqs.fastq
cat seq/splitRuns/SlOutRun1Euk/seqs.fna seq/splitRuns/SlOutRun2Euk/seqs.fna > seq/SlOutEuk/seqs.fna
cat seq/splitRuns/SlOutRun1Euk/seqs.fastq seq/splitRuns/SlOutRun2Euk/seqs.fastq > seq/SlOutEuk/seqs.fastq

#Make a UPARSE Directory
mkdir seq/UPARSEout
mkdir seq/Euk/UPARSEout

# get quality stats
usearch64 -fastq_stats seq/SlOut/seqs.fastq -log seq/UPARSEout/seqs.stats.log
usearch64 -fastq_stats seq/SlOutEuk/seqs.fastq -log seq/Euk/UPARSEout/seqs.stats.log

# remove low quality reads (trimming not required for paired-end data)
usearch64 -fastq_filter seq/SlOut/seqs.fastq -fastaout seq/UPARSEout/seqs.filtered.fasta -fastq_maxee 0.5
usearch64 -fastq_filter seq/SlOutEuk/seqs.fastq -fastaout seq/Euk/UPARSEout/seqs.filtered.fasta -fastq_maxee 5

# dereplicate seqs
usearch64 -derep_fulllength seq/UPARSEout/seqs.filtered.fasta  -fastaout seq/UPARSEout/seqs.filtered.derep.fasta -sizeout
usearch64 -derep_fulllength seq/Euk/UPARSEout/seqs.filtered.fasta  -fastaout seq/Euk/UPARSEout/seqs.filtered.derep.fasta -sizeout

# filter singletons
usearch64 -sortbysize seq/UPARSEout/seqs.filtered.derep.fasta -minsize 2 -fastaout seq/UPARSEout/seqs.filtered.derep.mc2.fasta
usearch64 -sortbysize seq/Euk/UPARSEout/seqs.filtered.derep.fasta -minsize 2 -fastaout seq/Euk/UPARSEout/seqs.filtered.derep.mc2.fasta

# cluster OTUs (de novo chimera checking can not be disabled in usearch)
usearch64 -cluster_otus seq/UPARSEout/seqs.filtered.derep.mc2.fasta -otus seq/UPARSEout/seqs.filtered.derep.mc2.repset.fasta
usearch64 -cluster_otus seq/Euk/UPARSEout/seqs.filtered.derep.mc2.fasta -otus seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.fasta

# reference chimera check
usearch64 -uchime_ref seq/UPARSEout/seqs.filtered.derep.mc2.repset.fasta -db ~/Fast_Analysis/QIIME/gold.fasta  -strand plus -nonchimeras seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta
usearch64 -uchime_ref seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.fasta -db ~/Fast_Analysis/QIIME/gold.fasta  -strand plus -nonchimeras seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta

# label OTUs using UPARSE python script
fasta_number.py seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta OTU_ > seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta
fasta_number.py seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta OTU_ > seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta

#Make an otus folder
mkdir otus/
mkdir otus/Euk/

#Copy this file to a repset.fna file for later use
cp seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta otus/RepSet.fna
cp seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta otus/Euk/RepSet.fna

# map the _original_ quality filtered reads back to OTUs
usearch64 -usearch_global seq/UPARSEout/seqs.filtered.fasta -db seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta -strand plus -id 0.97 -uc seq/UPARSEout/otu.map.uc
usearch64 -usearch_global seq/Euk/UPARSEout/seqs.filtered.fasta -db seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta -strand plus -id 0.97 -uc seq/Euk/UPARSEout/otu.map.uc

#Modify OTU table for input into QIIME
#HOLD UP
#I had to modify the uc2otutab.py script according to directions on the QIIME forums
#def GetSampleId(Label):
#	Fields = Label.split(";")
#	for Field in Fields:
#		if Field.startswith("barcodelabel="):
#			return Field[13:]
#	die.Die("barcodelabel= not found in read label '%s'" % Label)
#Present in the original script was modified to just 
#def GetSampleId(Label): 
#    SampleID = Label.split()[0].split('_')[0] 
#    return SampleID 
#Allowing the script to proceed.
python /home/lab/bin/uc2otutab.py seq/UPARSEout/otu.map.uc > seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt
python /home/lab/bin/uc2otutab.py seq/Euk/UPARSEout/otu.map.uc > seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt


# assign taxonomy 
echo “Assigning Taxonomy”
assign_taxonomy.py -m mothur -t /media/lab/Storage/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt -r /media/lab/Storage/Silva119_release/rep_set/97/Silva_119_rep_set97.fna -i otus/RepSet.fna -o otus/TaxonomyOut/
assign_taxonomy.py -m mothur -t /media/lab/Storage/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt -r /media/lab/Storage/Silva119_release/rep_set/97/Silva_119_rep_set97.fna -i otus/Euk/RepSet.fna -o otus/Euk/TaxonomyOut/

#SCP Files to a server that doesn't have issues with BIOM conversion. 
scp -i ~/gb2015.pem otus/TaxonomyOut/RepSet_tax_assignments.txt ubuntu@54.211.79.11:/home/ubuntu/Bac/
scp -i ~/gb2015.pem seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt ubuntu@54.211.79.11:/home/ubuntu/Bac/
scp -i ~/gb2015.pem otus/Euk/TaxonomyOut/RepSet_tax_assignments.txt ubuntu@54.211.79.11:/home/ubuntu/Euk/
scp -i ~/gb2015.pem seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt ubuntu@54.211.79.11:/home/ubuntu/Euk/
# convert to biom. Done on another server slightly differently. 
biom convert --table-type="OTU table" -i seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt -o otus/UPARSE.biom --to-hdf5
biom convert --table-type="OTU table" -i seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt -o otus/Euk/UPARSE.biom --to-hdf5

# add taxonomy to BIOM table. Done on the server with no BIOM issues. 
echo “Adding Metadata”
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/TaxonomyOut/RepSet_tax_assignments.txt -i otus/UPARSE.biom -o otus/otuTable.biom
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/Euk/TaxonomyOut/RepSet_tax_assignments.txt -i otus/Euk/UPARSE.biom -o otus/Euk/otuTable.biom

#Remove Euks from Bac analysis, and Bac from Euk Analysis
filter_taxa_from_otu_table.py -i Bac/otuTable.biom -o Bac/NoEuk.biom -n D_0__Eukaryota
filter_taxa_from_otu_table.py -i Euk/otuTable.biom -o Euk/NoBac.biom -n D_0__Bacteria

#Summarize BIOM files, and send back to the main server. These are decompressed, and placed in the correct folder.
biom summarize-table -i Bac/NoEuk.biom -o Bac/Summary.txt
biom summarize-table -i Euk/NoBac.biom -o Euk/Summary.txt
scp Bac.tar.gz lab@129.15.47.230:/media/lab/Work/Active_QIIME/B20/InitialAnalysis/
scp Euk.tar.gz lab@129.15.47.230:/media/lab/Work/Active_QIIME/B20/InitialAnalysis/

#Align and Filter Sequences
align_seqs.py -i otus/RepSet.fna -t /media/lab/Storage/Silva119_release_aligned_rep_files/97_16S_only/Silva_119_rep_set97_aligned_16S_only.fna -o otus/RepSet_Aligned/
filter_alignment.py -i otus/RepSet_Aligned/RepSet_aligned.fasta -o otus/RepSet_Aligned/ -e 0.10
make_phylogeny.py -i otus/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/RepSet.tre

align_seqs.py -i otus/Euk/RepSet.fna -t /media/lab/Storage/Silva119_release_aligned_rep_files/97_18S_only/Silva_119_rep_set97_aligned_18S_only.fna -o otus/Euk/RepSet_Aligned/
filter_alignment.py -i otus/Euk/RepSet_Aligned/RepSet_aligned.fasta -o otus/Euk/RepSet_Aligned/ -e 0.10
make_phylogeny.py -i otus/Euk/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/Euk/RepSet.tre


#Run Beta diversity analyses- seemingly what we're interested in at this point.
beta_diversity_through_plots.py -i otus/NoEuk.biom -m tags.txt -o BacBdiv500/ -t otus/RepSet.tre -e 500
beta_diversity_through_plots.py -i otus/Euk/NoBac.biom -m tags.txt -o EukBdiv300/ -t otus/Euk/RepSet.tre -e 300

