#!/bin/sh

# finds all animal names from the source directory 

find . -type f | sed -ne 's%.*trimmed/\(IT[[:digit:]]\{2\}\).*\.fna$%\1%pg' | sort -u

