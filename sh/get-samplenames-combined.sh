#!/bin/sh

# pushd () { command pushd "$@" > /dev/null; }
# popd () { command popd "$@" > /dev/null; }

PREFIX=/export/bse/reads/old

printf '# 454 samples\n'
find $PREFIX/454 -type f | sed -ne 's%.*trimmed/\(IT[[:digit:]]\{2\}\).*\.fna$%\1%pg' | sort -u
find $PREFIX/454 -type f | sed -ne 's%.*/\(RA[[:digit:]]\{2\}\).*\.fna$%\1%pg' | sort -u

printf '# illumina samples\n'
find $PREFIX/paul -type f | sed -ne 's%.*/\([[:digit:]]-Gr[[:digit:]]\)-[[:digit:]]\{2\}..*\.fastq$%\1%pg' | sort -u

printf '# amprolium samples\n'
find $PREFIX/illumina -type f | sed -ne 's%.*\([[:alpha:]]\{2\}[[:digit:]]\{2\}\).*\.fastq$%\1%pg' | sort -u
find $PREFIX/illumina -type f | sed -ne 's%.*\(S.-0-Gruppe.\).*\.fastq$%\1%pg' | sort -u

