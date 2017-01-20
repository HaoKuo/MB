#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
    usage: python duplex_pro.py ubam spacer output-name
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

#fastq_dir = "/work-z/user/guoh/tech-RD/molecular_barcode/fastq/"
#test_fq1 = "/work-z/user/guoh/tech-RD/molecular_barcode/fastq/S370_05B_CHG015512-Sample20161103-1-I13-AACGTGAT_L007_R1.fastq.gz"

#fastqToBam_dir = "/work-z/user/guoh/tech-RD/molecular_barcode/fastqTobam/"
#uBam_file = fastqToBam_dir + "S1_I13.unaligned.bam"


def consensus_caller(input_reads, cutoff, tag, length_check):

    nuc_identity_list = [0, 0, 0, 0, 0, 0]  # In the order of T, C, G, A, N, Total
    nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
    consensus_seq = ''
    if length_check is True:
        for read in input_reads[1:]:
            if len(read) != len(input_reads[0]):
                raise Exception("Read lengths for tag %s used for calculating the SSCS are not uniform!!!" % tag)
    for i in xrange(len(input_reads[0])):  # Count the types of nucleotides at a position in a read.
        for j in xrange(len(input_reads)):  # Do this for every read that comprises a tag family.
            try:
                if input_reads[j][i] == 'T':
                    nuc_identity_list[0] += 1
                elif input_reads[j][i] == 'C':
                    nuc_identity_list[1] += 1
                elif input_reads[j][i] == 'G':
                    nuc_identity_list[2] += 1
                elif input_reads[j][i] == 'A':
                    nuc_identity_list[3] += 1
                elif input_reads[j][i] == 'N':
                    nuc_identity_list[4] += 1
                else:
                    nuc_identity_list[4] += 1
                nuc_identity_list[5] += 1
            except:
                break
        try:
            for j in [0, 1, 2, 3, 4]:
                if float(nuc_identity_list[j])/float(nuc_identity_list[5]) >= cutoff:
                    consensus_seq += nuc_key_dict[j]
                    break
                elif j == 4:
                    consensus_seq += 'N'
        except:
            consensus_seq += 'N'
        nuc_identity_list = [0, 0, 0, 0, 0, 0]  # Reset for the next nucleotide position

    return consensus_seq

def qual_calc(qual_list):
    return [sum(qual_score) for qual_score in zip(*qual_list)]


@timeit
def main():

    #right_spacer = "AGTCAGCTA"
    #right_spacer = "ATCGACTGA"
    #right_spacer = "TAGCTGACT"
    #left_spacer = "TAGCTGACT"
    uBam_file = sys.argv[1]
    left_spacer = sys.argv[2]
    right_spacer = left_spacer
    output = sys.argv[3]
    in_bam = uBam_file
    in_bam_file = pysam.AlignmentFile(in_bam, "rb", check_sq=False)
    paired_end_count = 1
    l_counter = 0
    r_counter = 0

    l_counter1 = 0
    r_counter1 = 0

    l_counter2 = 0
    r_counter2 = 0

    l_pos_1=defaultdict(lambda: 0)
    r_pos_1=defaultdict(lambda: 0)
    l_pos_2 = defaultdict(lambda: 0)
    r_pos_2 = defaultdict(lambda: 0)
    counter1 = 0
    counter2 = 0

    i = 0
    for line in itertools.islice(in_bam_file.fetch(until_eof=True), 0, None):
        if line.query_alignment_sequence.find(left_spacer) + 1:
            l_counter += 1
        elif line.query_alignment_sequence.find(right_spacer) + 1:
            r_counter += 1
        #print line.query_name
        if paired_end_count % 2 == 1:
            counter1 += 1
            if line.query_alignment_sequence.find(left_spacer) + 1:
                l_counter1 += 1
                l_pos_1[line.query_alignment_sequence.find(left_spacer) + 1] +=1
            elif line.query_alignment_sequence.find(right_spacer) + 1:
                r_counter1 += 1
                r_pos_1[line.query_alignment_sequence.find(right_spacer) + 1] += 1

        if paired_end_count % 2 == 0:
            counter2 += 1
            if line.query_alignment_sequence.find(left_spacer) + 1:
                l_counter2 += 1
                l_pos_2[line.query_alignment_sequence.find(left_spacer) + 1] += 1
            elif line.query_alignment_sequence.find(right_spacer) + 1:
                r_counter2 += 1
                r_pos_2[line.query_alignment_sequence.find(right_spacer) + 1] += 1

        paired_end_count+=1
        i+=1

    f=open(output+"-spacer_count.txt","wb")
    f.write("Total reads in fastqToBam file : {}\n\
Total reads in fq1 file : {}\n\
Total reads in fq2 file : {}\n".format(i,counter1,counter2))
    f.write("Total left_spacer count: {}\n\
fq1 file left_spacer count: {}\n\
fq2 file left_spacer count: {}\n\
".format(l_counter,l_counter1,l_counter2))
    f.write("#####################################\n")
    for p in l_pos_1:
        f.write("left_spacer position in fq1 file: {}\tcount: {}\n".format(p,l_pos_1[p]))
    f.write("#####################################\n")
    for p in l_pos_2:
        f.write("left_spacer position in fq2 file: {}\tcount: {}\n".format(p,l_pos_2[p]))
    f.write("#####################################\n")
    #for p in r_pos_1:
    #    f.write("right_spacer position in fq1 file: {}\tcount: {}\n".format(p,r_pos_1[p]))
    #f.write("#####################################\n")
    #for p in r_pos_2:
    #    f.write("right_spacer position in fq2 file: {}\tcount: {}\n".format(p,r_pos_2[p]))

    f.close()
    #print l_pos_1, l_pos_2
        #print line.query_alignment_sequence.find("TAAA")
        #print line.query_alignment_qualities


@timeit
def foo():
    in_bam = uBam_file
    tag_len = 12
    spcr_len = 9
    tagstats = False
    minmem=3
    maxmem=200
    cutoff=0.7
    Ncutoff=1
    write_sscs=False
    without_dcs=False
    rep_filt = 9
    prefix = "S1_I13.dcs"

    dummy_header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1575, 'SN': 'chr1'}, {'LN': 1584, 'SN': 'chr2'}]}
    in_bam_file = pysam.AlignmentFile(in_bam, "rb", check_sq=False)
    temp_bam = pysam.AlignmentFile(prefix + ".temp.bam", 'wb', header=dummy_header)
    paired_end_count = 1
    if write_sscs is True:
        read1_sscs_fq_file = gzip.open(prefix + '_read1_sscs.fq.gz', 'wb')
        read2_sscs_fq_file = gzip.open(prefix + '_read2_sscs.fq.gz', 'wb')
    if without_dcs is False:
        read1_dcs_fq_file = gzip.open(prefix + '_read1_dcs.fq.gz', 'wb')
        read2_dcs_fq_file = gzip.open(prefix + '_read2_dcs.fq.gz', 'wb')
    for line in itertools.islice(in_bam_file.fetch(until_eof=True),0,1000000):
        if paired_end_count %2 ==1:
            temp_read1_entry = pysam.AlignedSegment()
            temp_read1_entry.query_name = line.query_name
            temp_read1_entry.query_sequence = line.query_alignment_sequence
            temp_read1_entry.query_qualities = line.query_alignment_qualities
        if paired_end_count % 2 == 0:
            temp_bam_entry = pysam.AlignedSegment()
            if temp_read1_entry.query_sequence[:tag_len] > line.query_alignment_sequence[:tag_len]:
                temp_bam_entry.query_name = temp_read1_entry.query_sequence[:tag_len] + line.query_alignment_sequence[:tag_len] + '#ab'

            elif temp_read1_entry.query_sequence[:tag_len] < line.query_alignment_sequence[:tag_len]:
                temp_bam_entry.query_name = line.query_alignment_sequence[:tag_len] + temp_read1_entry.query_sequence[:tag_len] + '#ba'

            elif temp_read1_entry.query_sequence[:tag_len] == line.query_alignment_sequence[:tag_len]:
                paired_end_count += 1
                continue
            # Write entries for Read 1
            temp_bam_entry.query_name += ":1"
            temp_bam_entry.query_sequence = temp_read1_entry.query_sequence[tag_len + spcr_len:]
            temp_bam_entry.query_qualities = temp_read1_entry.query_qualities[tag_len + spcr_len:]
            temp_bam_entry.set_tag('X?', temp_read1_entry.query_name, 'Z')
            temp_bam.write(temp_bam_entry)
            # Write entries for Read 2
            temp_bam_entry.query_name = temp_bam_entry.query_name.replace('1', '2')
            temp_bam_entry.query_sequence = line.query_sequence[tag_len + spcr_len:]
            temp_bam_entry.query_qualities = line.query_qualities[tag_len + spcr_len:]
            temp_bam_entry.set_tag('X?', line.query_name, 'Z')
            temp_bam.write(temp_bam_entry)

        paired_end_count += 1

    in_bam_file.close()
    temp_bam.close()
    print "Sorting reads on tag sequence..."
    pysam.sort("-n", prefix + ".temp.bam", "-o", prefix + ".temp.sort.bam")
    os.remove(prefix + ".temp.bam")
    seq_dict = {'ab:1': [], 'ab:2': [], 'ba:1': [], 'ba:2': []}
    qual_dict = {'ab:1': [], 'ab:2': [], 'ba:1': [], 'ba:2': []}
    fam_size_x_axis = []
    fam_size_y_axis = []

    read1_dcs_len = 0
    read2_dcs_len = 0
    in_bam_file = pysam.AlignmentFile(prefix + '.temp.sort.bam', "rb", check_sq=False)
    first_line = in_bam_file.next()
    seq_dict[first_line.query_name.split('#')[1]].append(first_line.query_sequence)
    qual_dict[first_line.query_name.split('#')[1]].append(list(first_line.query_qualities))
    tag_count_dict = defaultdict(lambda: 0)

    for line in in_bam_file.fetch(until_eof=True):
        tag,subtag_order= first_line.query_name.split('#')[0],first_line.query_name.split('#')[1]
        if line.query_name.split('#')[0] == tag:
            seq_dict[line.query_name.split('#')[1]].append(line.query_sequence)
            qual_dict[line.query_name.split('#')[1]].append(list(line.query_qualities))
        else:
            # if len(seq_dict['ab:1']) != len(seq_dict['ab:2']) or len(seq_dict['ba:1']) != len(seq_dict['ba:2']):
            #     raise Exception('ERROR: Read counts for Read1 and Read 2 do not match for tag %s' % tag)
            for tag_subtype in seq_dict.keys():
                if len(seq_dict[tag_subtype]) > 0:
                    tag_count_dict[len(seq_dict[tag_subtype])] += 1
                if len(seq_dict[tag_subtype]) < minmem:
                    seq_dict[tag_subtype] = []
                    qual_dict[tag_subtype] = []
                elif minmem<=len(seq_dict[tag_subtype])<=maxmem:

                    print seq_dict

                    seq_dict[tag_subtype] = [consensus_caller(seq_dict[tag_subtype], cutoff, tag, False),
                                             str(len(seq_dict[tag_subtype]))]
                    qual_dict[tag_subtype] = qual_calc(qual_dict[tag_subtype])

                elif len(seq_dict[tag_subtype]) > maxmem:
                    seq_dict[tag_subtype] = [consensus_caller(seq_dict[tag_subtype][:maxmem], cutoff, tag, False),
                                             str(len(seq_dict[tag_subtype]))]
                    qual_dict[tag_subtype] = qual_calc(qual_dict[tag_subtype])


            if write_sscs is True:

                if len(seq_dict['ab:1']) != 0 and len(seq_dict['ab:2']) != 0:
                    corrected_qual_score = map(lambda x: x if x < 41 else 41, qual_dict['ab:1'])
                    read1_sscs_fq_file.write('@%s#ab/1\n%s\n+%s\n%s\n' %
                                                 (tag, seq_dict['ab:1'][0], seq_dict['ab:1'][1], "".join(chr(x + 33) for x in corrected_qual_score)))

                    corrected_qual_score = map(lambda x: x if x < 41 else 41, qual_dict['ab:2'])
                    read2_sscs_fq_file.write('@%s#ab/2\n%s\n+%s\n%s\n' %
                                                 (tag, seq_dict['ab:2'][0], seq_dict['ab:2'][1], "".join(chr(x + 33) for x in corrected_qual_score)))

                if len(seq_dict['ba:1']) != 0 and len(seq_dict['ba:2']) != 0:
                    corrected_qual_score = map(lambda x: x if x < 41 else 41, qual_dict['ba:1'])
                    read1_sscs_fq_file.write('@%s#ba/1\n%s\n+%s\n%s\n' %
                                                 (tag, seq_dict['ba:1'][0], seq_dict['ba:1'][1], "".join(chr(x + 33) for x in corrected_qual_score)))

                    corrected_qual_score = map(lambda x: x if x < 41 else 41, qual_dict['ba:1'])
                    read2_sscs_fq_file.write('@%s#ba/2\n%s\n+%s\n%s\n' %
                                                 (tag, seq_dict['ba:2'][0], seq_dict['ba:2'][1], "".join(chr(x + 33) for x in corrected_qual_score)))
                        # i=0
    # handle = gzip.open(test_fq1, "r")
    # test_dic = SeqIO.index(handle,"fastq")
    # print len(test_dic)


   # for record in itertools.islice(SeqIO.parse(handle,"fastq"),0,10):
   #     i+=1
   #     print record
   # print i



# def main():
#     foo()





if __name__ == '__main__':
    main()
