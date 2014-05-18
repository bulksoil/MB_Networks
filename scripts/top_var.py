#!/usr/bin/python

"""
This script will take in an OTU table with the rows as the OTUs and the columns at samples
It will:
	1) Calculate the variances for each OTU
	2) Grab the top x OTUs with the highest variance, where x is a user defined value
	3) It will print these OTUs out
"""

import sys
import getopt
import re
import scipy.stats as stats

def main():
	in_file, out_file, keep = opt_load()
	otu_vars = otu_load(in_file)
	top_otus = var_sort(otu_vars, keep)
	otu_print(open(in_file), top_otus, out_file)

def opt_load():
	in_file = ''
	out_file = ''
	keep = 1000
	opts, args = getopt.getopt(sys.argv[1:], "i:o:k:")
	for o, a in opts:
		if o == '-i':
			in_file = a
		elif o == '-o':
			out_file = a
		elif o == '-k':
			keep = int(a)
		else:
			print script_info['usage']
			sys.exit(2)
	return(in_file, out_file, keep)


def otu_load(in_file):
	lines = open(in_file)
	header = lines.readline()
	otu_vars = {}
	for line in lines:
		counts = line.rstrip("\n").split("\t")
		otu = counts.pop(0)
		float_counts = [float(x) for x in counts]
		var = var_calc(float_counts)
		otu_vars[otu] = var
	return(otu_vars)

def var_calc(counts):
	var = 0
	mean = sum(counts) / len(counts)
	for val in counts:
		var += (mean - val) ** 2
	return(var)

def var_sort(otu_vars, keep):
	top_otus = {}
	iters = 0
	for otu in sorted(otu_vars, key = otu_vars.get, reverse = True):
		if iters == keep:
			break
		iters += 1
		top_otus[otu] = 1
	return(top_otus)

def otu_print(lines, top_otus, out_file):
	fout = open(out_file, 'w')
	print>>fout, lines.readline().rstrip("\n")
	for line in lines:
		otu = line.split("\t")[0]
		if otu in top_otus:
			print>>fout, line.rstrip("\n")

main()






