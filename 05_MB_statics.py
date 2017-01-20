#!/user/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
    usage: python MB_statics.py fqgz_dir/
'''
import os, sys, time, re, glob, itertools
import subprocess, multiprocessing, shlex
import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import pysam
from timeFunc import timeit
from argparse import ArgumentParser
from collections import defaultdict



@timeit
def main():
    file_dic = sys.argv[1]
    file_list = glob.glob(file_dic + "*.fq.gz")
    pool = multiprocessing.Pool(10)
    status = []
    for gf in file_list:
    #gf = "fastqToBam22/I7-Tag11_read1_dcs.fq.gz"
        result = pool.apply_async(count_N,(gf,))
        status.append(result)
    pool.close()
    pool.join()
    rc =0
    for st in status:
        if st.successful():
            rc +=0
        else:
            rc =+1
    print "status: %s"%(rc)


def count_reads_N():
    pass



def count_N(gf):
    count_f = open(gf+"-statics.txt","wb")
    i=0
    nuc = defaultdict(lambda: 0)
    with gzip.open(gf,"rb") as f:
        for line in itertools.islice(f,1,None,4):
            Nper = round(line.strip().count("N") / len(line.strip()),2)
            nuc[Nper] +=1
            i += 1
    tl = sorted(nuc.iteritems(), key=lambda d:d[0], reverse = False )
    count_f.write("#File: %s\n"%(gf))
    count_f.write("#Total reads count: %s\n"%(i))
    for k in tl:
        count_f.write("N percent:\t%s\tcount:\t%s\n"%(k[0],k[1]))



if __name__ == '__main__':
    main()
