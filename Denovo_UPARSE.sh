#!/bin/bash

#To analyze multiple Illumina MiSeq runs using the most recent primer set, a slightly different pipeline is required
#This includes analyzing the 16S and 18S amplified together independently. 
#Software/Databases Required
#QIIME 1.9.1
#USEARCH 8
#NOTE: This version of my workflow requires a 64 bit version of UPARSE due to the number of samples processed. 
#SILVA r123 formatted for use within QIIME. I followed the instructions at https://github.com/mikerobeson/Misc_Code/tree/master/SILVA_to_RDP. 


# Use QIIME to prepare your FASTQ files. No reverse primer is required for 18S reads. 
echo "Beginning, extracting barcodes and reads" 
extract_barcodes.py -f seq/Run1Seq/Run1_Joined.fastq  -a -m Run1.txt -l 12 -o seq/splitRuns/preppedRun1/
extract_barcodes.py -f seq/Run2Seq/Run2_Joined.fastq  -a -m Run2.txt -l 12 -o seq/splitRuns/preppedRun2/
extract_barcodes.py -f seq/Run3Seq/Run3_Joined.fastq  -a -m Run3.txt -l 12 -o seq/splitRuns/preppedRun3/

#Separate 
extract_barcodes.py -f seq/Run1Seq/Run1_Euk1.fastq  -l 12 -o seq/splitRuns/preppedRun1EukA/
extract_barcodes.py -f seq/Run1Seq/Run1_Euk2.fastq  -l 12 -o seq/splitRuns/preppedRun1EukB/
extract_barcodes.py -f seq/Run2Seq/Run2_Euk1.fastq  -l 12 -o seq/splitRuns/preppedRun2EukA/
extract_barcodes.py -f seq/Run2Seq/Run2_Euk2.fastq  -l 12 -o seq/splitRuns/preppedRun2EukB/
extract_barcodes.py -f seq/Run3Seq/Run3_Euk1.fastq  -l 12 -o seq/splitRuns/preppedRun3EukA/
extract_barcodes.py -f seq/Run3Seq/Run3_Euk2.fastq  -l 12 -o seq/splitRuns/preppedRun3EukB/

echo "Creating needed directories, and joining required files"
#I need to create separate prepped folders for the purpose of joining the R1 and R2 reads.
mkdir seq/splitRuns/preppedRun1Euk
mkdir seq/splitRuns/preppedRun2Euk
mkdir seq/splitRuns/preppedRun3Euk

#Now I need to make a SlOut folder for both the bacterial/archaeal analysis, and the eukaryotic one.
mkdir seq/SlOut/
mkdir seq/SlOutEuk/

#Then it's time to bring together the output of the eukaryotic extract_barcodes.py commands.
cat seq/splitRuns/preppedRun1EukA/reads.fastq seq/splitRuns/preppedRun1EukB/reads.fastq > seq/splitRuns/preppedRun1Euk/reads.fastq
cat seq/splitRuns/preppedRun1EukA/barcodes.fastq seq/splitRuns/preppedRun1EukB/barcodes.fastq > seq/splitRuns/preppedRun1Euk/barcodes.fastq
cat seq/splitRuns/preppedRun2EukA/reads.fastq seq/splitRuns/preppedRun2EukB/reads.fastq > seq/splitRuns/preppedRun2Euk/reads.fastq
cat seq/splitRuns/preppedRun2EukA/barcodes.fastq seq/splitRuns/preppedRun2EukB/barcodes.fastq > seq/splitRuns/preppedRun2Euk/barcodes.fastq
cat seq/splitRuns/preppedRun3EukA/reads.fastq seq/splitRuns/preppedRun3EukB/reads.fastq > seq/splitRuns/preppedRun3Euk/reads.fastq
cat seq/splitRuns/preppedRun3EukA/barcodes.fastq seq/splitRuns/preppedRun3EukB/barcodes.fastq > seq/splitRuns/preppedRun3Euk/barcodes.fastq

#Use QIIME to demultiplex the data, with -q 0. Store output as fastq format (we will quality filter with usearch). This part of the analysis was taken from suggestions within the QIIME forums on best practices to integrate UPARSE into a usable QIIME workflow.
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun1/reads.fastq -b seq/splitRuns/preppedRun1/barcodes.fastq -m Run1.txt --barcode_type 12 -o seq/splitRuns/SlOutRun1/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun2/reads.fastq -b seq/splitRuns/preppedRun2/barcodes.fastq -m Run2.txt --barcode_type 12 -o seq/splitRuns/SlOutRun2/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun3/reads.fastq -b seq/splitRuns/preppedRun3/barcodes.fastq -m Run3.txt --barcode_type 12 -o seq/splitRuns/SlOutRun3/

#Now for the eukaryotes.
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun1Euk/reads.fastq -b seq/splitRuns/preppedRun1Euk/barcodes.fastq -m Run1.txt --barcode_type 12 -o seq/splitRuns/SlOutRun1Euk/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun2Euk/reads.fastq -b seq/splitRuns/preppedRun2Euk/barcodes.fastq -m Run2.txt --barcode_type 12 -o seq/splitRuns/SlOutRun2Euk/
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/splitRuns/preppedRun3Euk/reads.fastq -b seq/splitRuns/preppedRun3Euk/barcodes.fastq -m Run3.txt --barcode_type 12 -o seq/splitRuns/SlOutRun3Euk/

#And then, I use cat to concatenate the output of the multiple split_libraries_fastq.py commands.
cat seq/splitRuns/SlOutRun1/seqs.fna seq/splitRuns/SlOutRun2/seqs.fna seq/splitRuns/SlOutRun3/seqs.fna > seq/SlOut/seqs.fna
cat seq/splitRuns/SlOutRun1/seqs.fastq seq/splitRuns/SlOutRun2/seqs.fastq seq/splitRuns/SlOutRun3/seqs.fastq > seq/SlOut/seqs.fastq
cat seq/splitRuns/SlOutRun1Euk/seqs.fna seq/splitRuns/SlOutRun2Euk/seqs.fna seq/splitRuns/SlOutRun3Euk/seqs.fna > seq/SlOutEuk/seqs.fna
cat seq/splitRuns/SlOutRun1Euk/seqs.fastq seq/splitRuns/SlOutRun2Euk/seqs.fastq seq/splitRuns/SlOutRun3Euk/seqs.fastq > seq/SlOutEuk/seqs.fastq

#Make a UPARSE Directory for each analysis
mkdir seq/UPARSEout
mkdir seq/Euk/
mkdir seq/Euk/UPARSEout

# get quality stats
usearch64 -fastq_stats seq/SlOut/seqs.fastq -log seq/UPARSEout/seqs.stats.log
usearch64 -fastq_stats seq/SlOutEuk/seqs.fastq -log seq/Euk/UPARSEout/seqs.stats.log

# remove low quality reads (trimming not required for paired-end data). Maxee is set high for the Eukaryotes due to the error prone region at the terminus of the unpaired read.
usearch64 -fastq_filter seq/SlOut/seqs.fastq -fastaout seq/UPARSEout/seqs.filtered.fasta -fastq_maxee 1 -threads 2
usearch64 -fastq_filter seq/SlOutEuk/seqs.fastq -fastaout seq/Euk/UPARSEout/seqs.filtered.fasta -fastq_maxee 5 -threads 2

# dereplicate seqs
usearch64 -derep_fulllength seq/UPARSEout/seqs.filtered.fasta  -fastaout seq/UPARSEout/seqs.filtered.derep.fasta -sizeout -threads 2
usearch64 -derep_fulllength seq/Euk/UPARSEout/seqs.filtered.fasta  -fastaout seq/Euk/UPARSEout/seqs.filtered.derep.fasta -sizeout -threads 2

# filter singletons
usearch64 -sortbysize seq/UPARSEout/seqs.filtered.derep.fasta -minsize 2 -fastaout seq/UPARSEout/seqs.filtered.derep.mc2.fasta -threads 2
usearch64 -sortbysize seq/Euk/UPARSEout/seqs.filtered.derep.fasta -minsize 2 -fastaout seq/Euk/UPARSEout/seqs.filtered.derep.mc2.fasta -threads 2

# cluster OTUs (de novo chimera checking can not be disabled in usearch)
usearch64 -cluster_otus seq/UPARSEout/seqs.filtered.derep.mc2.fasta -otus seq/UPARSEout/seqs.filtered.derep.mc2.repset.fasta
usearch64 -cluster_otus seq/Euk/UPARSEout/seqs.filtered.derep.mc2.fasta -otus seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.fasta

# reference chimera check
usearch64 -uchime_ref seq/UPARSEout/seqs.filtered.derep.mc2.repset.fasta -db ~/Fast_Analysis/QIIME/gold.fasta  -strand plus -nonchimeras seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta -threads 2
usearch64 -uchime_ref seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.fasta -db ~/Fast_Analysis/QIIME/gold.fasta  -strand plus -nonchimeras seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta -threads 2

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
usearch64 -usearch_global seq/UPARSEout/seqs.filtered.fasta -db seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta -strand plus -id 0.97 -uc seq/UPARSEout/otu.map.uc -threads 2
usearch64 -usearch_global seq/Euk/UPARSEout/seqs.filtered.fasta -db seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta -strand plus -id 0.97 -uc seq/Euk/UPARSEout/otu.map.uc -threads 2

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


# assign taxonomy using SILVA 123
echo “Assigning Taxonomy”
assign_taxonomy.py -m mothur -t /media/lab/Storage/Silva_123/Taxonomy.txt -r /media/lab/Storage/Silva_123/RepSet/97_RepSet.fasta -i otus/RepSet.fna -o otus/TaxonomyOut/
assign_taxonomy.py -m mothur -t /media/lab/Storage/Silva_123/Taxonomy.txt -r /media/lab/Storage/Silva_123/RepSet/97_RepSet.fasta -i otus/Euk/RepSet.fna -o otus/Euk/TaxonomyOut/

#SCP Files to a server that doesn't have issues with BIOM conversion. 
scp -i ~/gb2015.pem otus/TaxonomyOut/RepSet_tax_assignments.txt ubuntu@23.22.42.26:/home/ubuntu/Bac/
scp -i ~/gb2015.pem seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt ubuntu@23.22.42.26:/home/ubuntu/Bac/
scp -i ~/gb2015.pem otus/Euk/TaxonomyOut/RepSet_tax_assignments.txt ubuntu@23.22.42.26:/home/ubuntu/Euk/
scp -i ~/gb2015.pem seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt ubuntu@23.22.42.26:/home/ubuntu/Euk/
# convert to biom. Done on another server slightly differently. 
biom convert --table-type="OTU table" -i seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt -o otus/UPARSE.biom --to-hdf5
biom convert --table-type="OTU table" -i seq/Euk/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt -o otus/Euk/UPARSE.biom --to-hdf5

# add taxonomy to BIOM table. Done on the server with no BIOM issues. 
echo “Adding Metadata”
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/TaxonomyOut/RepSet_tax_assignments.txt -i otus/UPARSE.biom -o otus/otuTable.biom
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/Euk/TaxonomyOut/RepSet_tax_assignments.txt -i otus/Euk/UPARSE.biom -o otus/Euk/otuTable.biom

#Summarize BIOM files, and send back to the main server. These are decompressed, and placed in the correct folder.
biom summarize-table -i Bac/NoEuk.biom -o Bac/Summary.txt
biom summarize-table -i Euk/NoBac.biom -o Euk/Summary.txt
#Remove Euks from Bac analysis, and Bac from Euk Analysis
filter_taxa_from_otu_table.py -i otus/otuTable.biom -o otus/16S_Only.biom -n D_0__Eukaryota
filter_taxa_from_otu_table.py -i otus/Euk/otuTable.biom -o otus/18S_Only.biom -n D_0__Bacteria

#Align and Filter Sequences

align_seqs.py -i otus/Euk/RepSet.fna -t /media/lab/Storage/Silva_123/FilteredAlignments/Cleaned/97_RepSetAligned_Unfiltered_pfiltered_filtered.fasta -o otus/RepSet_Aligned/
filter_alignment.py -i otus/RepSet_Aligned/RepSet_aligned.fasta -o otus/RepSet_Aligned/ -e 0.10
make_phylogeny.py -i otus/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/RepSet.tre

align_seqs.py -i otus/Euk/RepSet.fna -t /media/lab/Storage/Silva_123/FilteredAlignments/Cleaned/97_RepSetAligned_Unfiltered_pfiltered_filtered.fasta -o otus/Euk/RepSet_Aligned/
filter_alignment.py -i otus/Euk/RepSet_Aligned/RepSet_aligned.fasta -o otus/Euk/RepSet_Aligned/ -e 0.10
make_phylogeny.py -i otus/Euk/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/Euk/RepSet.tre

#Remove samples not needed to run the specific analysis
#To simplify things, and remove some confounding samples that have no additional time points beyond a single "shiny rock"
#sampling event, I'm going to remove controls, single point samples, and old samples. This is for both 16S and 18S.

filter_samples_from_otu_table.py -i ../17Dec2015_Analysis/otus/18S_Only.biom -o OnlyField18S.biom --sample_id_fp samples_to_keep.txt
filter_samples_from_otu_table.py -i OnlyField18S.biom -o 18SAbove100.biom -n 100
filter_fasta.py -f ../17Dec2015_Analysis/otus/Euk/RepSet.fna -b otus/18SAbove100.biom -o 18S_Filtered.fna
filter_tree.py -i ../17Dec2015_Analysis/otus/Euk/RepSet.tre -f 18S_Filtered.fna -o 18S_Filtered.tre

filter_samples_from_otu_table.py -i ../17Dec2015_Analysis/otus/16S_Only.biom -o OnlyField16S.biom --sample_id_fp samples_to_keep.txt
filter_samples_from_otu_table.py -i OnlyField16S.biom -o 16SAbove100.biom -n 100
filter_fasta.py -f ../17Dec2015_Analysis/otus/RepSet.fna -b otus/16SAbove100.biom -o 16S_Filtered.fna
filter_tree.py -i ../17Dec2015_Analysis/otus/RepSet.tre -f 16S_Filtered.fna -o 16S_Filtered.tre

#Clean up intermediate files
rm Only*
mv *.tre otus/
mv *.fna otus/

#Run Beta diversity analyses- seemingly what we're interested in at this point. I will produce a stepwise beta diversity analysis at 100, 300, 500, 800, 1700, and 5000 sequences for 18S.
#16S will be 100, 500, 1000, 2500, 3500, and 7000.
beta_diversity_through_plots.py -i otus/16SAbove100.biom -m tags.txt -o 16S_Bdiv100/ -t otus/16S_Filtered.tre -e 100
beta_diversity_through_plots.py -i otus/16SAbove100.biom -m tags.txt -o 16S_Bdiv300/ -t otus/16S_Filtered.tre -e 500
beta_diversity_through_plots.py -i otus/16SAbove100.biom -m tags.txt -o 16S_Bdiv500/ -t otus/16S_Filtered.tre -e 1000
beta_diversity_through_plots.py -i otus/16SAbove100.biom -m tags.txt -o 16S_Bdiv800/ -t otus/16S_Filtered.tre -e 2500
beta_diversity_through_plots.py -i otus/16SAbove100.biom -m tags.txt -o 16S_Bdiv1700/ -t otus/16S_Filtered.tre -e 3500
beta_diversity_through_plots.py -i otus/16SAbove100.biom -m tags.txt -o 16S_Bdiv5000/ -t otus/16S_Filtered.tre -e 7000

beta_diversity_through_plots.py -i otus/18SAbove100.biom -m tags.txt -o 18S_Bdiv100/ -t otus/18S_Filtered.tre -e 100
beta_diversity_through_plots.py -i otus/18SAbove100.biom -m tags.txt -o 18S_Bdiv300/ -t otus/18S_Filtered.tre -e 300
beta_diversity_through_plots.py -i otus/18SAbove100.biom -m tags.txt -o 18S_Bdiv500/ -t otus/18S_Filtered.tre -e 500
beta_diversity_through_plots.py -i otus/18SAbove100.biom -m tags.txt -o 18S_Bdiv800/ -t otus/18S_Filtered.tre -e 800
beta_diversity_through_plots.py -i otus/18SAbove100.biom -m tags.txt -o 18S_Bdiv1700/ -t otus/18S_Filtered.tre -e 1700
beta_diversity_through_plots.py -i otus/18SAbove100.biom -m tags.txt -o 18S_Bdiv5000/ -t otus/18S_Filtered.tre -e 5000

#Compress these into something easily and quickly downloadable
tar cvzf BdivAnalyses_AFRL_28Dec2015B20.tar.gz 1*

echo "alpha_diversity:metrics equitability,goods_coverage,PD_whole_tree,observed_otus,kempton_taylor_q,enspie" > alpha_params.txt
alpha_rarefaction.py -p alpha_params.txt -n 50 -t otus/16S_Filtered.tre --min_rare_depth 100 -i otus/16SAbove100.biom -m tags.txt -o 16S_Arare/



