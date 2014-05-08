#!/usr/bin/python

import sys
import getopt
import re

def main():
	meta_file, out_file = opt_parse()
	names, genes, paths = meta_parse(open(meta_file))
	print_out(names, genes, paths, out_file)

def opt_parse():
	## Load the files in based on the command line arguments
	meta_file = ""
	out_file = ""
	opts, args = getopt.getopt(sys.argv[1:], "i:o:")
	for o, a in opts:
		if o == '-i':
			meta_file = a
		elif o == '-o':
			out_file = a
		else:
			print script_info['usage']
			sys.exit(2)
	return(meta_file, out_file)

def meta_parse(lines):
	names = name_parse(lines.readline())
	genes = gene_parse(names, lines.readline())
	paths = path_parse(names, lines.readline())
	return(names, genes, paths)


def name_parse(line):
	names = line.rstrip("\n").split("\t")
	names.pop()
	names.pop(0)
	return(names)

def gene_parse(names, line):
	genes = line.rstrip("\n").split("\t")
	genes.pop(0)
	gene_dict = {}
	for pos in range(0, len(names)):
		gene_dict[names[pos]] = re.sub('\[.*\]', '\t', genes[pos])

	return(gene_dict)

def path_parse(names, line):
	paths = line.rstrip("\n").split("\t")
	paths.pop(0)
	path_dict = {}
	for pos in range(0, len(names)):
		path = paths[pos]
		path = re.sub('\s+', '_', path)
		path_dict[names[pos]] = re.sub(';', '\t', path)

	return(path_dict)

def print_out(names, genes, paths, out_file):
	out = open(out_file, 'w')
	print>>out, "KO\tDescription\tPath1\tPath2\tPath3"
	for ko in names:
		gene = genes[ko]
		path = paths[ko]
		to_print = [ko, gene, path]
		print>>out, "\t".join(to_print)
	return()

main()