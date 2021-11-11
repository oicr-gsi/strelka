#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

echo ".vcf files:"
for f in *vcf.gz;do echo $f;zcat $f | cut -f 1-6 | md5sum;done
