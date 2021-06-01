#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find *vcf* -xtype f -size +0 | sed 's/.*\.//' | sort | uniq -c
