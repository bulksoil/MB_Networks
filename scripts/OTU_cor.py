#!/usr/bin/python

__author__ = "Joe Edwards"
__version__ = "1.0.0"
__maintainer__ = "Joe Edwards"
__email__ = "edwards@ucdavis.edu"


"""
This script will take in an OTU table which has several dependencies on format
	1) Row names must be OTU labels
	2) Must be tab delimited
	3) No quotes

Using the OTU table it will generate either pearson or spearman correlation values based on the method chosen
"""

import sys
import getopt 
import re
from scipy.stats.stats import *



def main():
	otu_counts_file, out_file, method, otu_by_row = opt_load()
	otu_counts, otus = otu_load(open(otu_counts_file), otu_by_row)
	otu_cor(otus, otu_counts, out_file, method)


def opt_load():
	## Load the files in based on the command line arguments
	otu_counts = ''
	out_file = ''
	method = 'pearson'
	otu_by_row = 1
	opts, args = getopt.getopt(sys.argv[1:], "i::o:m:c")
	for o, a in opts:
		if o == '-i':
			otu_counts = a
		elif o == '-o':
			out_file = a
		elif o == '-m':
			method = str(a)
		elif o == '-c':
			otu_by_row = 0
		else:
			print script_info['usage']
			sys.exit(2)
	return(otu_counts, out_file, method, otu_by_row)

def otu_load(lines, otu_by_row):
	if otu_by_row == 1:
		samples = lines.readline().rstrip("\n")
		otu_counts = {}
		otus = []
		for line in lines:
			counts = line.rstrip("\n").split("\t")
			otu = str(counts.pop(0))
			otus.append(otu)
			con_counts = [float(x) for x in counts]
			otu_counts[otu] = con_counts
		return(otu_counts, otus)

def otu_cor(otus, otu_counts, out_file, method):
	name_print = ["OTU"]
	for otu in otus:
		name_print.append(otu)
	print "\t".join(name_print)
	for otu1 in otus:
		to_print = [otu1]
		for otu2 in otus:
			if method == "pearson":
				R = pearsonr(otu_counts[otu1], otu_counts[otu2])[0]
			if method == "spearman":
				R = spearmanr(otu_counts[otu1], otu_counts[otu2])[0]
			to_print.append(R)

		print "\t".join([str(x) for x in to_print])

main()
