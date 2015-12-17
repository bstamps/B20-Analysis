#!/bin/bash

#Join Paired Ends
pear -f seq/R1.fastq -r seq/R2.fastq -o seq/Amp -p 0.001 -v 100 -m 450 -n 250 -y 500m -j 4

# Use QIIME to prepare your FASTQ files. 
extract_barcodes.py -f seq/Amp.assembled.fastq  -a -m tags.txt -l 12 -o seq/prepped/

# Use QIIME to demultiplex the data, with -q 0. Store output as fastq format (we will quality filter with usearch). This part of the analysis was taken from suggestions within the QIIME forums on best practices to integrate UPARSE into a usable QIIME workflow.
split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i seq/prepped/reads.fastq -b seq/prepped/barcodes.fastq -m tags.txt --barcode_type 12 -o seq/SlOut/

#Make a UPARSE Directory
mkdir seq/UPARSEout

# get quality stats
usearch64 -fastq_stats seq/SlOut/seqs.fastq -log seq/UPARSEout/seqs.stats.log

# remove low quality reads (trimming not required for paired-end data)
usearch64 -fastq_filter seq/SlOut/seqs.fastq -fastaout seq/UPARSEout/seqs.filtered.fasta -fastq_maxee 0.5

# dereplicate seqs
usearch64 -derep_fulllength seq/UPARSEout/seqs.filtered.fasta  -fastaout seq/UPARSEout/seqs.filtered.derep.fasta -sizeout

# filter singletons
usearch64 -sortbysize seq/UPARSEout/seqs.filtered.derep.fasta -minsize 2 -fastaout seq/UPARSEout/seqs.filtered.derep.mc2.fasta

# cluster OTUs (de novo chimera checking can not be disabled in usearch)
usearch64 -cluster_otus seq/UPARSEout/seqs.filtered.derep.mc2.fasta -otus seq/UPARSEout/seqs.filtered.derep.mc2.repset.fasta

# reference chimera check
usearch64 -uchime_ref seq/UPARSEout/seqs.filtered.derep.mc2.repset.fasta -db /home/ubuntu/SILVA119/gold.fasta -strand plus -nonchimeras seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta

# label OTUs using UPARSE python script
fasta_number.py seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.fasta OTU_ > seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta

#Make an otus folder
mkdir otus/

#Copy this file to a repset.fna file for later use
cp seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta otus/RepSet.fna

# map the _original_ quality filtered reads back to OTUs
usearch64 -usearch_global seq/UPARSEout/seqs.filtered.fasta -db seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTUs.fasta -strand plus -id 0.97 -uc seq/UPARSEout/otu.map.uc

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
python /usr/bin/uc2otutab.py seq/UPARSEout/otu.map.uc > seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt

# convert to biom
biom convert --table-type="OTU table" -i seq/UPARSEout/seqs.filtered.derep.mc2.repset.nochimeras.OTU-table.txt -o otus/UPARSE.biom --to-hdf5

# assign taxonomy 
echo “Assigning Taxonomy”
assign_taxonomy.py -m mothur -t /media/lab/Storage/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt -r /media/lab/Storage/Silva119_release/rep_set/97/Silva_119_rep_set97.fna -i otus/RepSet.fna -o otus/TaxonomyOut/

# add taxonomy to BIOM table
echo “Adding Metadata”
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/TaxonomyOut/RepSet_tax_assignments.txt -i otus/UPARSE.biom -o otus/otuTable.biom

#Sort the OTU table based on order desired of samples
sort_otu_table.py -i otus/otuTable.biom -o otus/sortedOTUtable.biom -m tags.txt -s SortOrder

#Computing Summaries of Taxa
echo "Computing Summaries"
biom summarize-table -i otus/sortedOTUtable.biom -o Summary.txt
biom summarize-table -i otus/sortedOTUtable.biom --qualitative -o QualSummary.txt

#Write out a Summary
summarize_taxa_through_plots.py -i otus/sortedOTUtable.biom -p ../qPCR_BAFB_34_UncoatedQIIME/plotSummaryprefs.txt -o TaxaSummary/
tar cvzf TaxaSummaryForCaitlynPoster1_18S_Only.tar.gz TaxaSummary/

#Align and Filter Sequences
align_seqs.py -i otus/RepSet.fna -t /media/lab/Storage/Silva119_release_aligned_rep_files/97_16S_only/Silva_119_rep_set97_aligned_16S_only.fna -o otus/RepSet_Aligned/
filter_alignment.py -i otus/RepSet_Aligned/RepSet_aligned.fasta -o otus/RepSet_Aligned/ -e 0.10
make_phylogeny.py -i otus/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/RepSet.tre


beta_diversity_through_plots.py -i otus/sortedOTUtable.biom -m tags.txt -o bdiv300/ -t otus/RepSet.tre -e 300
beta_diversity_through_plots.py -i otus/sortedOTUtable.biom -m tags.txt -o bdiv1k/ -t otus/RepSet.tre -e 1000
beta_diversity_through_plots.py -i otus/sortedOTUtable.biom -m tags.txt -o bdiv5k/ -t otus/RepSet.tre -e 5000
beta_diversity_through_plots.py -i otus/sortedOTUtable.biom -m tags.txt -o bdiv10k/ -t otus/RepSet.tre -e 10000

jackknifed_beta_diversity.py -i otus/sortedOTUtable.biom -o BdivJack300/ -e 300 -m tags.txt -t otus/RepSet.tre -a -O 2
jackknifed_beta_diversity.py -i otus/sortedOTUtable.biom -o BdivJack1K/ -e 1000 -m tags.txt -t otus/RepSet.tre -a -O 2
jackknifed_beta_diversity.py -i otus/sortedOTUtable.biom -o BdivJack5K/ -e 5000 -m tags.txt -t otus/RepSet.tre -a -O 2
jackknifed_beta_diversity.py -i otus/sortedOTUtable.biom -o BdivJack10K/ -e 10000 -m tags.txt -t otus/RepSet.tre -a -O 2

mkdir Bdiv
mv BdivJack* Bdiv/
mv bdiv* Bdiv/

tar cvzf BdivAFRL16S.tar.gz Bdiv/
#Done
echo “Completed! Happy QIIMEing”
