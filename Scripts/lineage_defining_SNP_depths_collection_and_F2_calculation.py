##########################################################################################################################################################
#Important Packages
import vcf
import os
import pandas as pd
import numpy as np
import ast
import itertools
import time
import sys
import pickle
import shutil
import subprocess as sb
from os import popen, system, listdir, getcwd,path

##########################################################################################################################################################



##########################################################################################################################################################

#INPUTS

#path to full-zipped VCF file
path_to_full_VCF = sys.argv[1]

#will overwrite directory specified if it already exists
output_dir = sys.argv[2]

#path to list of lineage-defining reference positions
lineage_def_ref_pos = sys.argv[3]

#isolate identifier
tag = sys.argv[4]

# PATH to ./Coll2014_LinSpeSNPs_final.csv
Coll2014_LinSpeSNPs_final_CSV_PATH = sys.argv[5]

##########################################################################################################################################################



##########################################################################################################################################################

#take Full VCF and create a reduced version that contains only the Reference Positions for Lineage-defining SNP sites

#create directory in F2 mixed calc. for Luca directory
#if os.path.exists(output_dir):
#	shutil.rmtree(output_dir)
#	os.makedirs(output_dir)
#elif not os.path.exists(output_dir):
#	os.makedirs(output_dir)

#change current working directory to one in which output will be stored in
#os.chdir(output_dir)

cmd1 = f"time bcftools view {path_to_full_VCF} -O b -o {output_dir}/current.bcf"
system(cmd1)

cmd2 = f"bcftools index {output_dir}/current.bcf"
system(cmd2)

cmd3 = f"bcftools view {output_dir}/current.bcf --regions-file {lineage_def_ref_pos} -o {output_dir}/{tag}_lineage_positions.vcf -O v"
system(cmd3)

reduced_VCF = f"{output_dir}/{tag}_lineage_positions.vcf"

##########################################################################################################################################################




##########################################################################################################################################################

#Lineage Defining SNPs - Import Lineage Defining SNP sets from Coll et. al. 2014
lineage_defining_SNPs = pd.read_csv(Coll2014_LinSpeSNPs_final_CSV_PATH)

#According to Wyllie et. al. 2018, drop SNP sets corresponding to branches with fewer than 20 SNPs
excluded_branches = ['lineage1.2' , 'lineage3.1' , 'lineage3.1.2' , 'lineage4.1.2' , 'lineage4.3.4.2.1' , 'lineage4.6' , 'lineage4.7']

#create a filter to filter out any SNP corresponding to the excluded lineages
included_branch_filter = [lineage not in excluded_branches for lineage in lineage_defining_SNPs.lineage]
lineage_defining_SNPs = lineage_defining_SNPs[included_branch_filter]
lineage_defining_SNPs.reset_index(inplace = True)

##########################################################################################################################################################




##########################################################################################################################################################
def check_call_is_SNP(ref_allele , alt_alleles):
	
	'''This function checks to see if Call is a SNP and not an InDel or Structural Variant'''
	
	#check that Reference Allele is 1 base
	if len(ref_allele) == 1:
		
		#check to see if there is no alternate allele
		if alt_alleles == [None]:
			
			good_SNP = True
			
		#if there is an alternate allele(s) check to see that they are all just 1 base
		elif ( sum( [(len(alt_allele) == 1) for alt_allele in alt_alleles] ) == len(alt_alleles) ):
			
			good_SNP = True
			
		#at least 1 alternate allele was longer than 1 base
		else:
			
			good_SNP = False
	
	#reference allele was longer than 1 base        
	else:
		good_SNP = False
			
	return good_SNP
##########################################################################################################################################################



##########################################################################################################################################################
def get_lineage_defining_SNP_depths(reduced_VCF , output_dir):

	'''
	This function takes as input, the path to the reduced VCF file, then extracts SINGLE base calls from
	the VCF file generated by Pilon specified as a lineage defining SNP from Coll et. al. 2014. 
	The function returns a DataFrame for a single isolate that contains base calling information
	for all of these SNP positions.
	'''

	#path to VCF file
	VCF_file = reduced_VCF

	vcf_reader = vcf.Reader(open(VCF_file , 'r'))

	#create dictionaries to store information for each call
	lineage_dict = {}
	position_dict = {}
	depth_dict = {}
	max_base_count_dict = {}

	#indexer for dataframe containing lineage defining SNPs
	index = 0

	#iterate through each record from VCF file
	for record in vcf_reader:

		#check to see that the record does not correspond to a variant call that is a Structural Variant (SVTYPE)
		if 'SVTYPE' not in record.INFO.keys():

			#position of variant on the Reference Genome
			reference_position = record.POS 

			#reference allele
			ref_allele = record.REF

			#alternate alleles (if there are any)
			alt_alleles = record.ALT

			#check to see if call is at a lineage defining site AND that call is a SNP
			if ( reference_position in list( lineage_defining_SNPs.position ) ) and check_call_is_SNP(ref_allele , alt_alleles):

				######## Retrieve Relevant information for filtering quality of Base Call ########
				# Mean Base Quality @ locus
				BQ = record.INFO['BQ']

				# Mean Mapping Quality @ locus
				MQ = record.INFO['MQ']

				# Number of Reads w/ Deletion 
				DC = record.INFO['DC']

				# Number of Reads w/ Insertion
				IC = record.INFO['IC']  
				
				# Depth of Valid Reads in Pileup
				VD = record.INFO['DP']

				### Filtering Criteria
				#---> Mean Base Quality > 20
				#---> Mean Mapping Quality > 30
				#---> No Reads Supporting Insertions
				#---> No Reads Supporting Deletions
				#---> Number of High Quality Reads >= 25
				if (BQ > 20) and (MQ > 30) and (DC == 0) and (IC == 0) and (VD >= 25): #SNP passed filtering criteria!

					# Valid read depth; some reads may have been filtered
					total_depth = record.INFO['DP']

					#Count of As, Cs, Gs, Ts at locus
					base_counts = record.INFO['BC']

					#most common base - depth
					most_common_base_depth = np.max(base_counts)

					#lineage that is defined by SNP
					lineage_defined = list( lineage_defining_SNPs[lineage_defining_SNPs.position == reference_position].lineage )[0]

					#After filtering for high-quality lineage-defining calls, store all of the pertinent information about the Base Call
					lineage_dict[index] = lineage_defined
					position_dict[index] = reference_position
					depth_dict[index] = total_depth
					max_base_count_dict[index] = most_common_base_depth

					index += 1

	#convert dictionaries to series
	lineage = pd.Series(lineage_dict)
	position = pd.Series(position_dict)
	depth = pd.Series(depth_dict)
	max_base_count = pd.Series(max_base_count_dict)

	#create DataFrame to hold all lineage defining SNPs
	lineage_SNP_depths_from_sample_DF = pd.DataFrame()
	lineage_SNP_depths_from_sample_DF['lineage'] = lineage
	lineage_SNP_depths_from_sample_DF['position'] = position
	lineage_SNP_depths_from_sample_DF['depth'] = depth
	lineage_SNP_depths_from_sample_DF['max_base_count'] = max_base_count

	#calculate the minor depth for each SNP site
	lineage_SNP_depths_from_sample_DF['minor_depth'] = lineage_SNP_depths_from_sample_DF['depth'] - lineage_SNP_depths_from_sample_DF['max_base_count']

	return lineage_SNP_depths_from_sample_DF
##########################################################################################################################################################



#Store the lineage depths in a pickled DF
##########################################################################################################################################################
#call functions to get lineage defining SNP depths for SNPs in all SNP sets
lineage_SNP_depths_from_sample_DF = get_lineage_defining_SNP_depths(reduced_VCF , output_dir)

#pickle DataFrame that has total depth & minor depth for all lineage defining SNPs for this isolate
lineage_SNP_depths_from_sample_DF.to_pickle(output_dir + '/' + tag + '_lineage_defining_SNP_depths.pkl')
##########################################################################################################################################################



#Calculate F2 measure from the depths at lineage defining sites
##########################################################################################################################################################
def calculate_minor_allele_fraction_per_SNP_set(lineage_SNPs_from_sample_DF):

	#get a list of all lineages (i.e. branches) that we have SNP information for
	all_lineages = list( set( list(lineage_SNPs_from_sample_DF.lineage) ) )

	#create a series that will store the minor allele fraction (p) for each lineage set
	minor_allele_fraction_per_SNP_set = pd.Series(index = all_lineages)

	for lineage in minor_allele_fraction_per_SNP_set.index:

		#subset to lineage-defining SNPs
		specific_lineage_SNPs_depths = lineage_SNPs_from_sample_DF[lineage_SNPs_from_sample_DF.lineage == lineage]

		#calculate the Total Depth across all SNP sites for this SNP set
		D = np.sum(specific_lineage_SNPs_depths.depth)

		#calculate the Total Minor Depth
		M = np.sum(specific_lineage_SNPs_depths.minor_depth)

		#calculate the minor allele fraction 
		p = float(M) / float(D)

		#store p for this SNP set in series
		minor_allele_fraction_per_SNP_set[lineage] = p

	#sort values in descending order
	minor_allele_fraction_per_SNP_set.sort_values(ascending = False , inplace = True)
	
	return minor_allele_fraction_per_SNP_set
##########################################################################################################################################################



##########################################################################################################################################################
def calculate_F2(sorted_minor_allele_fraction_per_SNP_set):

	#calculate F2 - get top 2 lineages with the highest minor allele frequencies
	F2_lineages = list( sorted_minor_allele_fraction_per_SNP_set[0:2].index )

	#subset to SNPs that belong to either SNP set in F2 lineages
	F2_lineage_filter = [lineage in F2_lineages for lineage in list(lineage_SNP_depths_from_sample_DF.lineage)]
	F2_lineage_SNP_depths_DF = lineage_SNP_depths_from_sample_DF[F2_lineage_filter]

	try:
		#calculate the Total Depth across all SNP sites for this SNPs belonging to F2 SNP sets
		D2 = np.sum(F2_lineage_SNP_depths_DF.depth)
		#calculate the Total Minor Depth
		M2 = np.sum(F2_lineage_SNP_depths_DF.minor_depth)
		#calculate the minor allele fraction
		p2 = float(M2) / float(D2)
	
	except AttributeError:
		# DF may be empty!
		p2 = np.nan
	
	return p2
##########################################################################################################################################################



#Store F2 measure for isolate in a txt file
##########################################################################################################################################################
#calculate minor allele frequency estimates for all SNP sets & sort in descending order
sorted_minor_allele_fraction_per_SNP_set = calculate_minor_allele_fraction_per_SNP_set(lineage_SNP_depths_from_sample_DF)

#calculate F2
F2 = calculate_F2(sorted_minor_allele_fraction_per_SNP_set)

file = open(output_dir + '/' + tag + '_F2.txt' , 'w')
file.write(str(F2))
file.close()
##########################################################################################################################################################