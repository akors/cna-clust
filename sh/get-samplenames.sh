#!/bin/sh

pushd () { command pushd "$@" > /dev/null; }
popd () { command popd "$@" > /dev/null; }

PREFIX=/export/bse/reads/old

printf '# 454 samples\n'
pushd $PREFIX/454
get-samplenames-from-454-dir.sh
popd

printf '# illumina samples\n'
pushd $PREFIX/paul
get-samplenames-from-illumina-dir.sh
popd

printf '# amprolium samples\n'
pushd $PREFIX/illumina
get-samplenames-from-amprolium-dir.sh
popd

