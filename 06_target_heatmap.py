#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    usage: python target_heatmap.py
'''
#from tools import *
#from tools.findVar import grep_gz_file, grep_chr_loc
import os, glob, sys, itertools, gzip
import pybedtools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import Series, DataFrame
from ggplot import *
import seaborn as sns
from timeFunc import timeit

from dplython import (DplyFrame, X, diamonds, select, sift, sample_n, \
     sample_frac, head, arrange, mutate, group_by, summarize, DelayFunction)
cwd = os.getcwd()
work_dir = cwd+"/"
mapped_dir = work_dir + "sorted_filtered/"

def main():
    detail_txt = "/work-z/shared/Mapping/sample20161220_25+18+18+25+25+22/mpileup/covdep_detail.txt.gz"
    data = pd.read_table(detail_txt,compression = 'gzip',sep='\t')
    sample_name = ["S1_I"+str(s) for s in xrange(1,7)]+["S2_I"+str(ss) for ss in xrange(26,32)]

    data1 = data[data['id'].isin(sample_name)]
    data2 = DplyFrame(data1) >> sift(X["cumPos"] !=0 ) >> select(X.id,X.cumPos,X.depth)
    data3 = data2 >> sift(X['depth']>10000)
    data3.to_csv("filtered_10k_panel6_pos_depth.tsv",index=False,sep='\t')
    g = sns.FacetGrid(data2, col="id", margin_titles=True, size=6,col_wrap=2)
    g.map(plt.scatter, "cumPos", "depth", color="#338844", s=10, edgecolor="white", lw=.25)
    g.set(xlim=(-1,1200000))
    g.savefig("MB_Sample"+"-Coverage_Uniformity.jpg")
    test()


@timeit
def test():
    map_info = pd.read_table("/work-z/user/guoh/tech-RD/MB4/sorted_filtered/all_flagstat.txt",sep = '\t')
    sample_name = ["S1_I"+str(s) for s in xrange(1,7)]+["S2_I"+str(ss) for ss in xrange(26,32)]
    map_info = map_info[map_info['sample'].isin(sample_name)]
    mapped_reads = map_info.set_index('sample')['mapped_reads']

    panel6 = pybedtools.BedTool('/work-z/user/guoh/Data/panel6.bed')
    bam_file_list = glob.glob(mapped_dir+"*I?.sorted.bam")+glob.glob(mapped_dir+"*I??.sorted.bam")
    bam_name_list = [os.path.basename(x).split(".")[0] for x in bam_file_list]
    header = ['Chr','Start','End','GeneSymbol','Transcript']+bam_name_list
    result_MultiCov = panel6.multi_bam_coverage(bams=bam_file_list).saveas(work_dir+"MultiCov.txt")
    panel6_multicov_file = work_dir+"MultiCov.txt"
    panel6_multicov_data = pd.read_table(panel6_multicov_file, sep = '\t', header = None)
    panel6_multicov_data.columns = header
    panel6_multicov_data.to_csv(work_dir+"panel6_multiCov.txt",index=False,sep='\t')
    panel6_multicov = panel6_multicov_data.iloc[:,5:]
    #sns.set_style("darkgrid", {'xtick.direction':90})
    #f, ax = plt.subplots(figsize=(16, 8))
    #g = sns.heatmap(panel6_multicov,cmap="Blues",yticklabels=False,xticklabels=True)
    #plt.xticks(rotation=90)
    #plt.savefig('panel6-heatmap.png')
    #plt.clf()
    p6Td = panel6_multicov.T.describe().T
    bed_file = panel6_multicov_data.iloc[:,:5]
    bed_describe = pd.concat([bed_file,p6Td],axis=1).sort(['mean','std'],ascending=[0,0])
    bed_describe.to_csv(work_dir+"panel6_multicov_describe.csv",index=False,sep='\t')
    panel6_multicov_nor = panel6_multicov.div(mapped_reads/1000000).dropna(axis=1,how='all')
    panel6_multicov_nor.round(decimals=2).to_csv("panel6_multicov_nor.tsv",index=False,sep='\t')
    sns.set(context="paper", style="white")
    f, ax = plt.subplots(figsize=(12, 9))
    cmap1 = sns.light_palette("#4CB391", as_cmap=True)
    sns.heatmap(panel6_multicov_nor,cmap=cmap1,vmax = 3000, yticklabels=False)
    plt.savefig('panel6-heatmap2.png')
    plt.clf()





if __name__ == "__main__":
    main()

