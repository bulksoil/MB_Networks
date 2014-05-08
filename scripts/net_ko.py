#!/usr/bin/python



__author__ = "Joe Edwards"
__version__ = "1.0.0"
__maintainer__ = "Joe Edwards"
__email__ = "edwards@ucdavis.edu"

import sys
import getopt 
import re
import scipy.stats as stats

"""
How this script works:
	1) Takes in a network file, needs to have the OTU as the first column and then the module as the last column
		- Read the network file, keep track of OTUs and which modules they belong in
	2) Removes OTUs that do not exist in the PICRUSt database
	3) Open KO description file. This has pathway information and gene information
	4) Read the file that links OTUs to KOs
		- Keep track of KO description counts in each module and as a total
	5) 
"""

"""
# To do
	Normalize for how many otus certain KOs came from within module
	In
"""

script_info = {}
script_info['brief'] = "This script will take in a network file generated from otu co-abundance data.  The first column needs to be the OTU and the last needs to be the moldule it belongs to"
script_info['usage'] = "usage: ./net_ko.py -n network.txt -k ko.txt"

# Read the network file to find good otus and modules
def main():
	## Load the files in based on the command line arguments
	network_file, ko_file, out_file, desc_file, exclude_mod = opt_load()
	print "---> Finding acceptable otus in network file"
	good_otus = good_otu_finder(open(network_file), exclude_mod)
	print "---> Looking up KO descriptions"
	ko_descs = ko_desc(open(desc_file), level = "paths")
	print "---> Parsing ko file"
	module_sums, whole_sums = ko_compile(open(ko_file), good_otus, ko_descs)
	print "---> Calculating hypergeometic probabilities for enrichment" 
	group_phyper(module_sums, whole_sums, ko_descs, out_file)


def opt_load():
	## Load the files in based on the command line arguments
	network_file = ''
	ko_file = ''
	out_file = ''
	desc_file = ''
	exclude = ""
	opts, args = getopt.getopt(sys.argv[1:], "n:k:o:d:x:")
	for o, a in opts:
		if o == '-n':
			network_file = a
		elif o == '-k':
			ko_file = a
		elif o == '-o':
			out_file = a
		elif o == '-x':
			exclude = str(a)
		elif o == '-d':
			desc_file = a
		else:
			print script_info['usage']
			sys.exit(2)
	return(network_file, ko_file, out_file, desc_file, exclude)

def good_otu_finder(lines, mod_exclude):
	header = lines.readline()
	useful_otus = {}
	discarded_otus = 0
	discarded_mod = 0
	for line in lines:
		line = line.rstrip("\n")
		fields = line.split("\t")
		otu = str(fields[0])
		module = str(fields.pop())
		if module == str(mod_exclude):
			discarded_mod += 1
			continue
		elif re.match('.*New.*', otu, re.I):
			discarded_otus += 1
			continue
		else:
			otu = re.sub('Otu', '', otu)
			otu = re.sub('\"', '', otu)
			useful_otus[str(otu)] = str(module)
	print "	Discarded", discarded_otus, "OTUs."
	print "	Discarded ", discarded_mod, "OTUs from Module", mod_exclude
	print "	Kept", len(useful_otus), "OTUs."


	return(useful_otus)

def ko_compile(ko_lines, otu_module, ko_descs):
	found = 0;
	processed = 0
	whole_counts = {}
	module_counts = {}
	to_go = len(otu_module)
	ko_list = ko_lines.readline().rstrip("\n").split("\t")
	ko_list.pop(0)
	ko_list.pop()

	for line in ko_lines:
		counts = line.rstrip("\n").split("\t")
		otu = str(counts.pop(0))
		if (processed % 10000) == 0:
			print "Processed", processed,":	Found", found
		if found == to_go:
			print "	All good OTUs found."
			break
		elif otu in otu_module:
			found += 1
			processed += 1
			counts.pop()
			module = str(otu_module[otu])
			if module not in module_counts:
				module_counts[module] = {}
			for val in range(0, len(ko_list)):
				function = ko_descs[ko_list[val]]
				if function in module_counts[module]:
					module_counts[module][function] += float(counts[val])
				else:
					module_counts[module][function] = float(counts[val])
				if function in whole_counts:
					whole_counts[function] += float(counts[val])
				else:
					whole_counts[function] = float(counts[val])
		else:
			processed += 1
			continue

		

	return(module_counts, whole_counts)






def ko_desc(ko_desc_file, level = "paths"):
	descriptions = {}
	for line in ko_desc_file:
		fields = line.rstrip("\n").split("\t")
		ko = fields.pop(0)
		desc = ""
		if level == "paths":
			desc = fields.pop()
		descriptions[ko] = desc
	return(descriptions)

def ko_totaller(totals, names, values):
	"""
		Keeps track of the total counts for each KO
	"""

	for x in range(0, len(values)):
				ko = names[x]
				totals[ko] += float(values[x])
	return totals

def module_ko_adder(module, current_sums, names, values):
	"""
		Keeps track of the KO sums in each module
	"""

	if module in current_sums:
		for x in range(0, len(names)):
			current_sums[module][names[x]] += float(values[x])
	else:
		current_sums[module] = {}
		for x in range(0, len(names)):
			current_sums[module][names[x]] = float(values[x])
	return(current_sums)

def compact_desc_counts(module_total, whole_total, desc):
	module_comp = {}
	whole_comp = {}
	for module in module_total:
		current = module_total[module]
		module_comp[module] = {}
		for ko in current:
			description = desc[ko]
			if description in module_comp[module]:
				module_comp[module][description] += int(current[ko])
			else:
				module_comp[module][description] = int(current[ko])

	for ko in desc:
		if ko in whole_total:
			description = desc[ko]
			if description in whole_comp:
				whole_comp[description] += whole_total[ko]
			else:
				whole_comp[description] = whole_total[ko]

	return(module_comp, whole_comp)


def group_phyper(module_total, whole_total, desc, out_file):
	"""
		This will take the numbers from the KO counts and do
		the p value crunching at different specified levels
	"""
	fout = open(out_file, 'w')
	print>>fout, "Module\tCategory\twhole_sum\tmodule_sum\twhole_cat\tmod_cat\tpvalue"
	to_go = len(module_total)
	done = 0
	for module in module_total:
		done += 1
		if module in module_total:
			current = module_total[module]
			print "	Working on module", module, done,"/",to_go
			mod_sum = 0
			for descr in current:
				mod_sum += int(current[descr])
			for cat in current:
				if current[cat] == 0:
					continue
				else:
					whole_sum = sum(whole_total.values())
					whole_cat = whole_total[cat]
					module_cat = current[cat]
					p_val = stats.hypergeom.sf(module_cat, whole_sum, whole_cat, mod_sum)
					print>>fout, module, "\t", cat, "\t", whole_sum, "\t", mod_sum, "\t", whole_cat, "\t", module_cat, "\t", p_val
		else:
			print "Can't find module", module
				
def p_adjust(list, method = "BH"):
	n = len(list)




main()












