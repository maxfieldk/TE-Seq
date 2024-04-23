#!/bin/bash
indir=$1
find $indir -name "*.bam" -exec sh -c 'modkit pileup "$0" "$0".bedmethyl.bed --ref $1/$(echo "$0" | awk -F "." "{print \$4\".cons.ref.fa\"}") --combine-strands --cpg' {} $indir  \;
