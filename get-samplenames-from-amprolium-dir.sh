#!/bin/sh

find . -type f | sed -ne 's%.*/\([[:digit:]]-Gr[[:digit:]]\)-[[:digit:]]\{2\}..*\.fastq$%\1%pg' | sort -u

