#!/usr/bin/env python
import csv
import pysam
import textwrap
import sys
import argparse

parser = argparse.ArgumentParser(description="Given a BED file of repeat regions, extract the appropriate sequence from a fasta file. The output is written to the screen, so the user is advised to use redirection.")
parser.add_argument("BED", help="BED file containing region coordinates.")
parser.add_argument("Fasta", help="Fasta file")
args = parser.parse_args()

if(args.BED == None or args.Fasta == None) :
    parser.print_help()
    sys.exit()

bed = csv.reader(open(args.BED, "r"), dialect="excel-tab")
fasta = pysam.Fastafile(args.Fasta)
of = sys.stdout
for line in bed :
    seq = fasta.fetch(line[0], int(line[1]), int(line[2]))
    of.write(">%s\n" % line[3])
    for l in textwrap.wrap(seq, 80) :
        of.write("%s\n" % l)
fasta.close()
