#!/bin/bash


#ls I{4..9}* | while read x; do  
#	a=${x#*/}
#	b=${a%%.*}
#	echo 'python UnifiedConsensusMaker.py --input '$x' --taglen 12 --spacerlen 9 --write-sscs --prefix '$b 
#	echo 'python UnifiedConsensusMaker.py --input '$x' --taglen 11 --spacerlen 9 --write-sscs --prefix '$b'-Tag11' 
#	echo 'python MB_classifier.py --input '$x' --taglen 12 --spacerlen 9 --write-sscs --prefix '$b 
#done | parallel 

#ls I10* | while read x; do  
#	a=${x#*/}
#	b=${a%%.*}
#	echo 'python UnifiedConsensusMaker.py --input '$x' --taglen 12 --spacerlen 9 --tagstats --write-sscs --prefix '$b 
	#echo 'python UnifiedConsensusMaker.py --input '$x' --taglen 11 --spacerlen 9 --write-sscs --prefix '$b'-Tag11' 
#	echo 'python MB_classifier.py --input '$x' --taglen 12 --spacerlen 9 --write-sscs --prefix '$b 
#done | parallel 


ls *unaligned.bam | while read x; do  
	a=${x#*/}
	b=${a%%.*}
	echo 'python Refined_UCM.py --input '$x' --taglen 12 --spacerlen 9 --tagstats --write-sscs --prefix '$b 
	#echo 'python Refined_UCM.py --input '$x' --taglen 11 --spacerlen 9 --tagstats --write-sscs --prefix '$b'-Tag11' 
#	echo 'python MB_classifier.py --input '$x' --taglen 12 --spacerlen 9 --write-sscs --prefix '$b 
done | parallel -j 12

