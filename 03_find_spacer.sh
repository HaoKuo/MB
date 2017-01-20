#!/bin/bash
ls fastqToBam/*.bam | while read x; do
	a=${x#*/}
	b=${a%%.*}
	echo 'python duplex_pro.py '$x' AGTCAGTCA spacer_dist_right/'$b
done | parallel -j 12
#python duplex_pro.py fastqToBam/I4.unaligned.bam TAGCTGACT spacer_dist/I4
