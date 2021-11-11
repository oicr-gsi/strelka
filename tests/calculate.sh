#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

echo ".vcf files:"
find *vcf* -xtype f -size +0 | grep -v tbi | xargs md5sum | sort -V
