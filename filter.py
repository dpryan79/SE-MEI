#!/usr/bin/env python
import csv
import argparse
import sys

parser = argparse.ArgumentParser(description="Take a sorted list of clusters and filter them such that the output contains only cluster bounds of clusters with more than one member.")
parser.add_argument("ifile", metavar="Input", help="Input BED file")
parser.add_argument("type", metavar="Repeat type", help="Repeat type label")
parser.add_argument("ofile", metavar="Output", help="Output BED file")
args = parser.parse_args()

if(args.ifile == None or args.ofile == None or args.type==None) :
    parser.print_help()
    sys.exit()

of = open(args.ofile, "w")
lastID = None
lastChr = None
lastStart = None
lastStop = None
lastCount = 0
for line in csv.reader(open(args.ifile, "r"), dialect="excel-tab") :
    if(line[3] != args.type) :
        continue

    if(line[4] == lastID) :
        lastCount += 1
        lastStop = line[2]
    else :
        if(lastCount > 1) :
            of.write("%s\t%s\t%s\n" % (lastChr, lastStart, lastStop))
        lastChr = line[0]
        lastStart = line[1]
        lastStop = line[2]
        lastID = line[4]
        lastCount = 1
if(lastCount > 1) :
    of.write("%s\t%s\t%s\n" % (lastChr, lastStart, lastStop))
of.close()
