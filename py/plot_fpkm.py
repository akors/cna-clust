#!/usr/local/bin/python3

import csv
import matplotlib.pyplot
import numpy

with open('genes.fpkm_tracking', 'rt') as infile:
    infile.readline() # skip header
    fpkm = [float(x[9]) for x in csv.reader(infile, delimiter='\t')]



(n, b, p) = matplotlib.pyplot.hist(fpkm, bins=20, log=True)
print("Bins:")
print(b)
print("Counts:")
print(n)
print("Total: ", len(fpkm))
print("Max: ", numpy.max(fpkm))

matplotlib.pyplot.show()

