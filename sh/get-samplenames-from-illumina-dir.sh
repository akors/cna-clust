#!/bin/sh

find . -type f | sed -ne 's%.*\([[:alpha:]]\{2\}[[:digit:]]\{2\}\).*\.fastq$%\1%pg' | sort -u
find . -type f | sed -ne 's%.*\(S.-0-Gruppe.\).*\.fastq$%\1%pg' | sort -u

