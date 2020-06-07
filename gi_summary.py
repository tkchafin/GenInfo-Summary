#!/usr/bin/python

import re
import sys
import os
import getopt
import sys
import time
import pandas as pd
import numpy as np
from Bio import Entrez
from collections import Counter

def main():
	params = parseArgs()

	#read accessions from input table
	hits = readAccessionTable(params.input, params.col)
	hits=hits[int(params.headers):] #remove header line
	print(hits)
	
	#filter missing values and track how many there were
	badchars = ["-"]
	before=len(hits)
	hits=[value for value in hits if value != "-"]
	tot=len(hits)
	nbad=before-tot
	
	print("Found",str(nbad), "lines missing an accession.")

	#get counts for each hit
	hits = sorted(hits)
	hit_count = Counter(hits)
	hits=hit_count.keys()
	uniq_count=len(hits)
	print("Found",str(tot),"total accessions, with",uniq_count,"being unique")

	#initialize table
	print("Fetching taxon names from accessions...")
	sci_names=efetchNamesFromAccessions(hits, params.batchSize, params.email, params.time, params.api)
	sci_names = sorted(sci_names)
	
	#get counts for each species name
	sci_counts=Counter(sci_names)
	print("Found",str(len(sci_counts.keys())),"unique taxa represented.")
	
	#NCBI query to convert tax names to taxids
	print("Fetching taxids...")
	taxids=efetchTaxidsFromNames(sci_counts.keys(), params.email, params.api)
	
	#Taxonomy query to get full taxonomies for each species
	print("Fetching complete taxonomy from taxids...")
	full_taxonomy=efetchTaxonomyFromTaxid(taxids.values(), params.email, params.api)
	
	#final counts
	#make a dataframe maybe, make sure to keep counts from sci_counts

def efetchTaxonomyFromTaxid(taxids, email, api=None):
	ret=dict()
	wait=0.33
	if api:
		wait=0.1
	
	#Specify email for Entrex session
	Entrez.email = email
	
	#post NCBI query
	for t in taxids:
		search_handle     = Entrez.efetch("taxonomy", id=t, retmode="xml")
		search_results    = Entrez.read(search_handle)
		lineage = {d['Rank']:d['ScientificName'] for d in search_results[0]['LineageEx']}
		print(lineage)
		time.sleep(wait)
		sys.exit()
	return(taxids)

#returns a dict mapping scientific (species) name to taxids
def efetchTaxidsFromNames(names, email, api=None):
	taxids = dict()
	
	wait=0.33
	if api:
		wait=0.1
	
	#Specify email for Entrex session
	Entrez.email = email

	#post NCBI query
	for t in names:
		search_handle     = Entrez.esearch("taxonomy", term=t)
		search_results    = Entrez.read(search_handle)
		taxids[t] = search_results['IdList'][0]
		time.sleep(wait)
	print(taxids)
	return(taxids)


def efetchNamesFromAccessions(hits, batchSize, email, waitTime, api=None):
	
	sci_names=list()
	
	#Specify email for Entrex session
	Entrez.email = email

	#Report how many hits were found
	if (len(hits) <= 0):
		print("No accessions found!\n")
		sys.exit(0)
	print("Found %s entries! Downloading in batches of %s...\n"%(len(hits), batchSize))

	#post NCBI query
	#print(hits)
	search_handle     = Entrez.epost("nucleotide", id=",".join(hits))
	search_results    = Entrez.read(search_handle)
	webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]

	#fetch batches of results, parse each batch
	count=1

	for start in range(0,len(hits),batchSize):
		print("Parsing batch %s..." %count)
		try:
			my_args = {'db':"nucleotide", 'retmode':'xml', 'retstart':start, 'retmax':batchSize, 'webenv':webenv, 'query_key':query_key}
			if api:
				my_args['api']=api
			handle = Entrez.efetch(**my_args)
			records = Entrez.parse(handle)
	
		except Exception as e:
			print("Oh no! Unknown exception on batch",count,":",e)
			count += 1
		else:
			#If no exception, parse records
			for rec in records:
				#print("rec")
				#print(rec["GBSeq_organism"].replace(' ', '+').strip())
				name=rec["GBSeq_organism"].split(" ")
				nj=name[0] + "+" + name[1]
				sci_names.append(nj)
				#print(len(rec["GBSeq_taxonomy"].split("; ")),rec["GBSeq_taxonomy"])
		count += 1
		#Wait for pre-defined time before attempting next query
		time.sleep(waitTime)
	return(sci_names)	

def readAccessionTable(tab, col):
	if col > 0:
		col = col-1
	else:
		col=1
		print("Warning: readAccessionTable(): Column number should not be less than 1. Setting to 1.")
	
	with open(tab, "r") as acc_lines:
		return([line.strip().split("\t")[col] for line in acc_lines])

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'i:e:A:t:b:hc:', \
			["email=","retmax=","batch=","out=", "input=",
			"help", "time=", "api=", "col=", "headers="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.email=None
		self.api=None
		self.input=None
		self.retmax=1000000000
		self.batchSize=500
		self.out=""
		self.col=1
		self.time = 0.5
		self.headers=1

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg in options:
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 'e' or opt == 'email':
				self.email = arg
			elif opt == 'i' or opt=='input':
				self.input=arg
			elif opt == 'h' or opt == 'help':
				pass
			elif opt == 'o' or opt == 'out':
				self.out = arg
			elif opt == 'm' or opt == 'retmax':
				self.retmax = int(arg)
			elif opt == 'b' or opt == 'batch':
				self.batchSize = int(arg)
			elif opt == 't' or opt == 'time':
				self.time = float(arg)
			elif opt == "A" or opt=="api":
				self.api=arg
			elif opt=="headers":
				self.headers=int(arg)
			elif opt=="c" or opt=="col":
				self.col=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.email:
			self.display_help("Email not provided <-e,--email>")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ngi_summary.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Summarizes taxonomy representation among sets of GIs/Accession numbers")
		print("""
	Mandatory Arguments:
		-i, --input	: [Required] Input table
		-e, --email	: [Required] Email used for Entrez
		
	Optional arguments:
		-c,--col	: Column number (1=first) containing GIs in input (default=1)
		-A, --api	: API key (needed for >3 queries per second to avoid HTTP timeout)
		-t,--time	: Time (seconds) to wait between batches (default=0.5)
		-b,--batch	: Batch size for Entrez queries [default=500]
		--headers	: Number of header lines to skip (default=1)
		-h,--help	: Displays help menu
""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
