#!/usr/bin/python
#!python

#Module setup
import sys, getopt
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Trim FASTQ files and reconfigure nexus-specific barcodes')
	parser.add_argument('-i', '--infile', required=True, help='Input FASTQ file')
	parser.add_argument('-o', '--outfile', required=True, help='Output FASTQ file')
	parser.add_argument('-t', '--trim', default=0, help='Trim all reads to a specific length before processing.')
	parser.add_argument('-k', '--keep', default=10, help='Minimum number of genomic bases needed to keep fragment.')
	parser.add_argument('-b', '--barcode', default='CTGA', help='Fixed barcode sequence(s) (comma-separated, no spaces) to filter for.')
	parser.add_argument('-r', '--randomidlength', default=5, help='Length of random barcode sequences that act as identifiers in later QC.')
	args = parser.parse_args()

	print(args)
	def process_fastq_chunk(fastq_chunk, args):
		print('pending')