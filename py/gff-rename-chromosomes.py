#!/usr/local/bin/python3

import sys
import csv

def load_reverse_nametable(filename):
    with open(filename, 'rt') as infile:
        table = {row[1] : row[0] for row in csv.reader(infile, delimiter=' ')}

    return table

def rewrite_table(tablename, infilename, outfilename):
    table = load_reverse_nametable(tablename)

    with open(infilename, 'rt') as infile, open(outfilename, 'wt') as outfile:
        for line in infile:

            if len(line) == 0 or line[0] == "#":
                outfile.write(line)
                continue

            fields = line.split("\t")
            fields[0] = table[fields[0]]
            outfile.write("\t".join(fields))



if __name__ == "__main__":
    rewrite_table(sys.argv[1], sys.argv[2], sys.argv[3])

