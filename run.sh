#!/bin/bash


ls I{4..9}* | while read x; do  
	a=${x#*/}
	b=${a%%.*}
	echo 'python UnifiedConsensusMaker.py --input '$x' --taglen 12 --spacerlen 9 --write-sscs --prefix '$b 
done | parallel 
