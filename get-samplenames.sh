#!/bin/sh

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

printf '# 454 samples\n'
pushd /export/bse/454
get-samplenames-from-454-dir.sh
popd

printf '# illumina samples\n'
pushd /export/bse/paul
get-samplenames-from-illumina-dir.sh
popd

printf '# amprolium samples\n'
pushd /export/bse/illumina
get-samplenames-from-amprolium-dir.sh
popd

