#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
    usage: python run_pipeline.py
'''
import os, sys, time, re, glob, itertools
import subprocess, multiprocessing, shlex
import gzip
from findVar import grep_chr_loc, grep_gz_file
from Bio.Seq import Seq
from Bio import SeqIO
import pysam
from guess_seq import guess_barcode

work_dir = os.getcwd()+"/" 

fastq_dir = "/work-z/data/projects/sample20161220_25+18+18+25+25+22/"
fastq_info = fastq_dir+"fastq_files.txt.gz"

tmp_dir = work_dir+"tmp/"
mapped_dir = work_dir + "mapped2/"
fastQC_dir = work_dir + "fastqc/"
guesseq_dir = work_dir + "guesseq/"
fastqToBam_dir = work_dir + "fastqToBam/"
fragment_dir = work_dir + "fragment/"
mpileup_all_dir = work_dir + "mpileup_all/"
pyvacheck_dir = work_dir + "pyvacheck/"
varscan_dir = work_dir + "varscan/"


for d in [tmp_dir, mapped_dir, fastQC_dir, fastqToBam_dir, guesseq_dir, mpileup_all_dir, varscan_dir, pyvacheck_dir, fragment_dir]:
    if not os.path.exists(d):
        os.mkdir(d)
    else:
        continue

VARSCAN = "/work-a/app/varscan-2.4.2/VarScan.v2.4.2.jar"
GATK = "/work-z/user/guoh/soft/GenomeAnalysisTK.jar"
PICARD = "/work-z/user/guoh/soft/picard/build/libs/picard.jar"
SAMTOOLS = "/work-z/user/guoh/soft/samtools-1.3.1/samtools"
BWA = "/work-a/app/bwa-0.7.12/bwa"
hg19_fasta = "/work-a/public/gatk_bundle/2.8/hg19/ucsc.hg19.fasta"


def main():

    pool = multiprocessing.Pool(12)
    
    sample_name = ["S1_I"+str(s) for s in xrange(1,7)]+\
            ["S2_I"+str(ss) for ss in xrange(26,32)]
    status = []
    for name in sample_name:
        lane = name.split("_")[0]
        i=name.split("_")[1]
        tmpl = grep_chr_loc(lane,i,fastq_info)
        L = tmpl[0].strip().split("\t")
        print lane, i 
        R1 = fastq_dir+L[3]
        R2 = fastq_dir+L[4]
        cmd1 = call_BwaMem(R1,R2,lane,i)
        cmd2 = call_FastQC(fastQC_dir, R1)
        cmd3 = call_FastQC(fastQC_dir, R2)
        cmd4 = call_FastqToSam(R1,R2,fastqToBam_dir+i+".unaligned.bam",i)
        CMD_list = [cmd1,cmd2,cmd3,cmd4]
        result = pool.apply_async(check_call_CMD, (CMD_list,))
        status.append(result)    
    pool.close()
    pool.join()
    rc = 0
    for st in status:
        if st.successful():
            rc += 0
        else:
            rc += 1
    if rc >0:
        print("Pipeline Failed! ")
    else:
        print("Done! Running Successfully!")
        status=[]
        pool = multiprocessing.Pool(12)
        bamFiles = glob.glob(mapped_dir+"*.bam")
        for bam in bamFiles:
            result=pool.apply_async(guess_barcode, args=(bam,guesseq_dir+os.path.basename(bam),1000000))
            status.append(result)
        pool.close()
        pool.join()
        pool = multiprocessing.Pool(12)
        f= glob.glob(mapped_dir + "*.bam")
        pool.map(run_one_cmd_pipeline, f)
        pool.close()
        pool.join()

def check_call_CMD(CMD_list):
    rc_l = []
    for cmd in CMD_list:
        command_line = ["/bin/bash", "-c", cmd]
        rc =  subprocess.check_call(command_line)
        rc_l.append(rc)
    return rc_l


def run_one_cmd_pipeline(bamFileAbsPath):
    filename = os.path.basename(bamFileAbsPath)
    sampleName = filename.split(".")[0]
    #print sampleName
    os.chdir(mapped_dir)
    cmd1 = call_SortBam(bamFileAbsPath, sampleName)
    cmd2 = call_BamFilter(sampleName + ".sorted.bam", sampleName)
    cmd3 = call_BamFlagstat(sampleName + ".filtered.sorted.bam")
    cmd4 = call_mpileup_all(sampleName + ".filtered.sorted.bam", mpileup_all_dir + sampleName)
    cmd5 = call_pyvarcheck(mpileup_all_dir + sampleName + ".mpileup.gz", pyvacheck_dir + sampleName + ".pyvarcheck.gz")
    cmd6 = call_varscan(mpileup_all_dir + sampleName + ".mpileup.gz", varscan_dir + sampleName, "snp")
    cmd7 = call_varscan(mpileup_all_dir + sampleName + ".mpileup.gz", varscan_dir + sampleName, "indel")
    CMD_list = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7]
    rc_list = []

    for i, pcmd in enumerate(CMD_list):
        command_line = ["/bin/bash", "-c", pcmd]
        #print ' '.join(command_line)
        rc = subprocess.check_call(command_line)
        #rc = popen_call(command_line)
        rc_list.append(rc)
    return rc_list


def run_workflow(mapped_dir):
    # bamFile="S1_I1.bam"
    # sampleName="S1_I1"
    # print call_SortBam(bamFile,sampleName)
    # print call_BamFilter(bamFile,sampleName)
    # print call_BamFlagstat(bamFile)
    # print call_mpileup_all(bamFile,sampleName)
    # print call_pyvarcheck("S1_I1.mpileup.gz","S1_I1.pyvarcheck.gz")
    os.chdir(mapped_dir)
    pool = multiprocessing.Pool(5)
    status = []
    for f in glob.glob(mapped_dir + "*.bam"):
        result = pool.apply_async(run_one_cmd_pipeline, (f,))
        status.append(result)
    pool.close()
    pool.join()
    rc = 0
    for st in status:
        if st.successful():
            rc += 0
        else:
            rc += 1
    if rc >0:
        print("Pipeline Failed! ")
    else:
        print("Done! Running Successfully!")
    return rc

    ##########
    ##  call varscan mpileup2snp and mpilup2indel
    ##########
    # rc_snp = call_varscan_mp(barcode_mpileup_datadic,"snp")
    # #log.write("return code of snp calling is {}\n".format(rc_snp))
    # if rc_snp == 0:
    #     log.write( "snp calling successfully\n")
    #     rc_indel = call_varscan_mp(barcode_mpileup_datadic, "indel")
    #     if rc_indel == 0:
    #         log.write("indel calling successfully\n")
    #     else: log.write("indel calling failed\n")
    # else: log.write("snp calling failed\n")
    #
    ##########
    ##  call varscan somatic
    ##########
    #
    # rc_somatic = call_varscan_somatic_mp(somatic_pair_dic_tbc)
    # log.write("return code of somatic is {}\n".format(rc_somatic))
    # if rc_somatic == 0: log.write("somatic calling successfully\n")
    # else: log.write("somatic calling failed\n")
    #
    ##########


def call_CombinVariants_mp(pppd1_snpindel_dic,option):
    pool = multiprocessing.Pool(9)
    status = []
    for c in pppd1_snpindel_dic:
        tmpdic = pppd1_snpindel_dic[c]
        tmpCMD = call_GATK_CombineVariants(tmpdic,somatic_merge_dir+c+".somatic.vcf",option)
        #print tmpCMD
        command_line = ['/bin/bash', '-c', tmpCMD]
        result = pool.apply_async(subprocess.check_call, (command_line,))
        status.append(result)
    pool.close()
    pool.join()
    rc =0
    for st in status:
        if st.successful():
            rc+=0
        else:
            rc+=1
    return rc


def popen_call(command_line):
    args = shlex.split(command_line)
    p = subprocess.Popen(args,stdout=None)
    stdinfo = p.communicate()
    return stdinfo


def call_varscan_mp(barcode_mpileup_datadic, type):
    pool = multiprocessing.Pool(9)
    status = []
    for b in barcode_mpileup_datadic:
        tmpPath = barcode_mpileup_datadic[b]
        tmpCMD = call_varscan(tmpPath, b, type)
        command_line = ['/bin/bash', '-c', tmpCMD]
        result = pool.apply_async(subprocess.check_call, (command_line,))
        status.append(result)
    pool.close()
    pool.join()
    rc =0
    for st in status:
        if st.successful():
            rc+=0
        else:
            rc+=1
    return rc


def call_varscan_somatic_mp(somatic_pair_dic_tbc):
    os.chdir(somatic_dir)
    pool = multiprocessing.Pool(9)
    status = []
    for b in somatic_pair_dic_tbc:
        tmpPathList = somatic_pair_dic_tbc[b]
        Normal_gz = tmpPathList[1]
        Tumor_gz = tmpPathList[0]
        tmpCMD = call_varscan_somatic(Normal_gz,Tumor_gz,b)
        command_line = ['/bin/bash', '-c', tmpCMD]
        result = pool.apply_async(subprocess.check_call, (command_line,))
        status.append(result)
    pool.close()
    pool.join()
    rc =0
    for st in status:
        if st.successful():
            rc+=0
        else:
            rc+=1
    os.chdir(work_dir)
    return rc


def call_varscan(mpileup_gz,fileBaseName, mtype):
    cmd1 = "java -Xmx32g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=" + tmp_dir + " -jar " + VARSCAN + " mpileup2" + mtype + " "
    cmd2 = "<(zcat " + mpileup_gz + ") "
    cmd3 = "--min-coverage 100 --min-reads2 4 --min-var-freq 0.001 --p-value 0.01 --strand-filter 1 --output-vcf --variants "
    cmd4 = ">" + fileBaseName + "." + mtype + ".vcf"
    cmd = cmd1 + cmd2 + cmd3 + cmd4
    return cmd


def call_varscan_somatic(normal_mpileup_gz, tumor_mpileup_gz, fileBaseName):
    cmd1 = "java -Xmx8g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=" + tmp_dir + " -jar " + varscan + " somatic "
    cmd2 = "<(zcat " + normal_mpileup_gz + ") "
    cmd3 = "<(zcat " + tumor_mpileup_gz + ") "
    cmd4 = fileBaseName
    ## tumor and normal sample coverage >10 reads
    cmd5 = " --min-coverage 10 --strand-filter 1 --min-var-freq 0.001 --somatic-p-value 0.05 --output-vcf"
    ## coverage >100 reads
    #cmd5 = " --min-coverage 10 --min-coverage-normal 100 --min-coverage-tumor 100 --strand-filter 1 --somatic-p-value 0.05 --output-vcf"
    cmd = cmd1 + cmd2 + cmd3 + cmd4 + cmd5
    return cmd


def call_GATK_CombineVariants(vcf_dic, output_path, merge_option):
    cmd1 = "java -Xmx8g -jar " + GATK + " -T CombineVariants "
    cmd2 = " -R " + hg19_fasta + " "
    tmp = []
    for k in vcf_dic:
        string = "--variant:%s %s"%(k,vcf_dic[k])
        tmp.append(string)
    cmd3 = " ".join(tmp)
    cmd4 = " -o " + output_path + " "
    cmd5 = "-genotypeMergeOptions " + merge_option
    cmd = cmd1 + cmd2 + cmd3 + cmd4 + cmd5
    return cmd

def call_SortBam_mp(bam_dic):
    pool = multiprocessing.Pool(9)
    status = []
    for b in bam_dic:
        tmpPath =bam_dic[b]
        tmpCMD = call_SortBam(bam_dic[b], b)
        command_line = ['/bin/bash', '-c', tmpCMD]
        result = pool.apply_async(subprocess.check_call, (command_line,))
        status.append(result)
    pool.close()
    pool.join()
    rc =0
    for st in status:
        if st.successful():
            rc+=0
        else:
            rc+=1
    return rc


def call_SortBam(bamFile, sampleName):
    cmd1 = "java -Xmx32g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=" + tmp_dir + " -jar " + PICARD + " SortSam "
    cmd2 = "INPUT=" + bamFile
    cmd3 = " OUTPUT=" + sampleName + ".sorted.bam CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=" + tmp_dir
    return cmd1 + cmd2 + cmd3


def call_BamFilter(bamFile, sampleName):
    cmd = SAMTOOLS +" view -F 0x900 " + bamFile + " -b -o " + sampleName + ".filtered.sorted.bam"
    return cmd


def call_BamFlagstat(bamFile):
    cmd = SAMTOOLS + " flagstat " + bamFile + " > " + bamFile + ".flagstat.txt"
    return cmd

def call_mpileup_all(bamFile, sampleName):
    cmd = SAMTOOLS + " mpileup -B -f " + hg19_fasta + " -q 1 -d 10000000 " + bamFile + " | gzip -9 > " + sampleName  + ".mpileup.gz"
    return cmd

def call_pyvarcheck(mpileup_gz,output_gz):
    cmd = "pypy /work-a/user/guoh/tools/pyVarcheck.py -i " + mpileup_gz + " -o " + output_gz
    return cmd


def call_FastqToBam(F1,F2,OUTPUT,SM):
    cmd1 = "java -Xmx32g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=" + tmp_dir + " -jar " + PICARD + " FastqToSam "
    cmd2 = "F1=" + F1 + " F2=" + F2 +" O="+OUTPUT+" SM="+SM
    return cmd1+cmd2

def call_BwaMem(R1,R2,lane,name):
    cmd1 = BWA + " mem -t 12 -M "+hg19_fasta+ " "+R1+" "+R2 +" -R "
    cmd2 = r"'@RG\tID:"+lane+"_"+name+r"\tSM:"+lane+"_"+name+r"\tLB:"+lane+"_"+name+r"\tPL:ILLUMINA\tPU:flowcell-barcode.lane'"
    cmd3 = "|"+SAMTOOLS+" view -Sb - > mapped2/"+lane+"_"+name+".bam"
    return cmd1+cmd2+cmd3

def call_FastQC(fastQC_dir,fastq_file):
    cmd1 = "fastqc -t 6 -o "+fastQC_dir+" "+fastq_file
    return cmd1

def call_FastqToSam(R1, R2, output, SM):
    cmd1 = "java -Xmx32g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=" + tmp_dir + " -jar "+PICARD+" FastqToSam " + "FASTQ="+R1+" FASTQ2="+R2+" OUTPUT="+output+" SM="+SM 
    return cmd1

if __name__ == '__main__':
    os.chdir(work_dir)
    main()

