##################################################################################################################################
#                                                database_registration_processors
# 
# primer_database_parser
# gene_fasta_generator
# amplicon_seq_generator
# exon_fasta_generator
# coding_seq_generator
#
##################################################################################################################################

import os, sys, glob, re, string
import numpy as np
import pandas as pd
from difflib import SequenceMatcher
from collections import defaultdict
from functions import seq_fast_generator, overlap_seq_finder

##################################################################################################################################






########################################################primer_database_parser###################################################

def primer_database_parser (refGene_path, master_file2in, master_file2out):

    genome_list = ['hg38', 'mm10', 'macFas5']
    gene2info = {}
    for genome_build in genome_list:
        refGene_file = genome_build + '_refGene.txt'
        with open(os.path.join(refGene_path, refGene_file), 'r') as ref2in:
            for line in ref2in:
                line = line.rstrip('\n')
                eles = line.split('\t')
                gene_id = genome_build + '-' + eles[1]
                gene_start = eles[2] + '_' + eles[3] + '_' + eles[4]
                gene_end = eles[2] + '_' + eles[3] + '_' + eles[5]
                exon_starts = eles[9]
                exon_ends = eles[10]
                gene_name = eles[12]
                gene_info = gene_start + '\t' + gene_end + '\t' + exon_starts + '\t' + exon_ends + '\t' + gene_name
                gene2info[gene_id] = gene_info

    dir_path = os.path.dirname(os.path.realpath(__file__))
    newlines = []
    with open(os.path.join(dir_path, master_file2in), 'r') as master2in:
        header = master2in.readline()
        header = header.rstrip('\n')
        for line in master2in:
            line = line.rstrip('\n')
            eles = line.split('\t')
            gene_id = eles[1] + '-' + eles[2]
            gene_info = gene2info[gene_id]
            newline = line + '\t' + gene_info + '\n'
            newlines.append(newline)

    with open(os.path.join(dir_path, master_file2out), 'w') as master2out:
        master2out.write(header + '\tgene_start\tgene_end\texon_starts\texon_ends\tgene_name_copy\n')
        for newline in newlines:
            master2out.write(newline)

#################################################################################################################################





            
########################################################gene_fasta_generator#####################################################

def gene_fasta_generator (master_file2in):

    dir_path = os.path.dirname(os.path.realpath(__file__))
    gene2seq = {}
    uniq = []

    with open(os.path.join(dir_path, master_file2in), 'r') as master2in:
        header = master2in.readline()
        for line in master2in:
            line = line.rstrip('\n')
            eles = line.split('\t')
            genome_build = eles[1]
            gene_id = eles[2]
            gene_name = eles[3]
            genome_gene = genome_build + '_' + gene_id + '_' + gene_name
            cor_start = eles[6]
            cor_end = eles[7]
            cor_start_eles = cor_start.split('_')
            cor_end_eles = cor_end.split('_')
            chrom = cor_start_eles[0]
            strand = cor_start_eles[1]
            gene_start = int(cor_start_eles[2])
            gene_end = int(cor_end_eles[2])

            path_fa = genome_build + '/fa_by_chrom/'
            path_fa_full = dir_path.replace('scripts', path_fa)

            if genome_gene not in uniq:
                seq = ''
                if 'chr' in chrom and '_' not in chrom:
                    if strand == '+':
                        seq = seq_fast_generator(path_fa_full, chrom, gene_start, gene_end)
                        seq = seq.upper()
                    elif strand == '-':
                        seq = seq_fast_generator(path_fa_full, chrom, gene_start, gene_end)
                        seq = seq.upper()
                        seq = seq[::-1]
                        seq = seq.translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))
                gene2seq[genome_gene] = seq
                uniq.append(genome_gene)

    with open(os.path.join(dir_path, 'gene_target_list.fa'), 'w') as fa2out:
        for genome_gene, seq in gene2seq.items():
            fa2out.write('>' + genome_gene + '\n' + seq + '\n')

#################################################################################################################################






########################################################amplicon_seq_generator###################################################

def amplicon_seq_generator (master_file2in, master_file2out):

    dir_path = os.path.dirname(os.path.realpath(__file__))

    with open(os.path.join(dir_path, 'gene_target_list.fa'), 'r') as genefa2in:
        geneid2fa = {}
        for line in genefa2in:
            line = line.rstrip('\n')
            if '>' in line:
                gene_id = line.replace('>', '')
            else:
                geneid2fa[gene_id] = line
    
    newlines = []
    with open(os.path.join(dir_path, master_file2in), 'r') as master2in:
        header = master2in.readline()
        header = header.rstrip('\n')
        for line in master2in:
            line = line.rstrip('\n')
            eles = line.split('\t')
            gid = eles[0]
            gene_fa_id = eles[1] + '_' + eles[2] + '_' + eles[3]
            gene_fa = geneid2fa[gene_fa_id]
            primer_f = eles[4]
            primer_r = eles[5]
            primer_frt = eles[4][::-1]
            primer_rrt = eles[5][::-1]
            primer_frt = primer_frt.translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))
            primer_rrt = primer_rrt.translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))

            mark_f = mark_r = ''
            pre_forward = post_forward = pre_reverse = post_reverse = amp_seq = ''
            if primer_f in gene_fa:
                mark_f = '+'
                pre_forward, post_forward = gene_fa.split(primer_f, 1)
                if primer_rrt in post_forward:
                    mark_r = '-'
                    pre_reverse, post_reverse = post_forward.split(primer_rrt, 1)
                    
                    #regex1 = re.compile(re.escape(primer_f) + r'.*' + re.escape(primer_rrt))
                    #regex_list1 = re.findall(primer_rrt, post_forward)
                    #print(regex_list1)
                    
                    amp_seq = primer_f + pre_reverse + primer_rrt
                else:
                    mark_r = 'primer_not_found'
            elif primer_frt in gene_fa:
                mark_f = '-'
                pre_forward, post_forward = gene_fa.split(primer_frt, 1)
                
                if primer_r in pre_forward:
                    mark_r = '+'
                    pre_reverse, post_reverse = pre_forward.split(primer_r, 1)
                    
                    regex_list2 = re.findall(primer_r, pre_forward)
                    
                    amp_seq = primer_r + post_reverse + primer_frt
                else:
                    mark_r = 'primer_not_found'
            else:
                mark_f = 'primer_not_found'
            
            ###########hard_code_for_cyno_ttr#########################
            if 'G0003'  in gid or 'G0004' in gid:
                amp_seq = 'GTCACTCCTACCTCATTTAGCGTGCATTTTAAATGTAGGAGCGGGATGTCACAGAAACACTCACCGTAGGACCAGCCTCAGACACAAATACCAGTCCAGCGAGGCAGAGGAGGAGCAGACGATGAGAAGCCATCCTGCCAAGAACGAGTGGACTTCTGTGATGGCTGCTC'
            if 'G0005'  in gid:
                amp_seq = 'GGGAGCAGCCATCACAGAAGTCCACTCATTCTTGGCAGGATGGCTTCTCATCGTCTGCTCCTCCTCTGCCTTGCTGGACTGGTATTTGTGTCTGAGGCTGGCCCTACGGTGAGTGTTTCTGTGACATCCCATTCCTACATTTAAGATTCACGCTAAATGAAGTAGAAGTGACTCCTTCCAGCTTTGCCAACCAGCTTTTA'
            ##########################################################

            newline = line + '\t' + amp_seq + '\t' + mark_f + '\t' + mark_r + '\n'
            newlines.append(newline)

    with open(os.path.join(dir_path, master_file2out), 'w') as master2out:
        master2out.write(header + '\tamplicon_sequence\tmark_forward_primer\tmark_reverse_primer\n')
        for newline in newlines:
            master2out.write(newline)

#################################################################################################################################


 



########################################################exon_fasta_generator#####################################################

def exon_fasta_generator (master_file2in):

    dir_path = os.path.dirname(os.path.realpath(__file__))
    exon2seq = {}
    uniq = []
    exon_starts = []
    exon_ends = []
    with open(os.path.join(dir_path, master_file2in), 'r') as master2in:
        header = master2in.readline()
        for line in master2in:
            line = line.rstrip('\n')
            eles = line.split('\t')
            genome_build = eles[1]
            gene_id = eles[2]
            gene_name = eles[3]
            genome_gene = genome_build + '_' + gene_id + '_' + gene_name
            cor_start = eles[6]
            cor_start_eles = cor_start.split('_')
            chrom = cor_start_eles[0]
            strand = cor_start_eles[1]
            path_fa = genome_build + '/fa_by_chrom/'
            path_fa_full = dir_path.replace('scripts', path_fa)

            exon_starts = eles[8].split(',')
            exon_ends = eles[9].split(',')
            exon_count = len(exon_starts) - 1
            if genome_gene not in uniq:
                for i in range(0, exon_count):
                    exon_id = genome_gene + '_exon_' + exon_starts[i] + '_' + exon_ends[i]
                    exon_start = int(exon_starts[i])
                    exon_end = int(exon_ends[i])
                    seq = ''
                    if 'chr' in chrom and '_' not in chrom:
                        if strand == '+':
                            seq = seq_fast_generator(path_fa_full, chrom, exon_start, exon_end)
                            seq = seq.upper()
                        elif strand == '-':
                            seq = seq_fast_generator(path_fa_full, chrom, exon_start, exon_end)
                            seq = seq.upper()
                            seq = seq[::-1]
                            seq = seq.translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))
                    exon2seq[exon_id] = seq
                uniq.append(genome_gene)

    with open(os.path.join(dir_path, 'gene_target_exons_list.fa'), 'w') as fa2out:
        for exon_id, seq in exon2seq.items():
            fa2out.write('>' + exon_id + '\n' + seq + '\n')

#################################################################################################################################






########################################################coding_seq_generator#####################################################

def coding_seq_generator (master_file2in, master_file2out):

    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'gene_target_exons_list.fa'), 'r') as exonfa2in:
        exonid2fa = {}
        for line in exonfa2in:
            line = line.rstrip('\n')
            if '>' in line:
                exon_id = line.replace('>', '')
            else:
                exonid2fa[exon_id] = line
    
    newlines = []
    with open(os.path.join(dir_path, master_file2in), 'r') as master2in:
        header = master2in.readline()
        header = header.rstrip('\n')
        for line in master2in:
            line = line.rstrip('\n')
            eles = line.split('\t')
            gid = eles[0]
            genome_build = eles[1]
            gene_id = eles[2]
            gene_name = eles[3]
            genome_gene = genome_build + '_' + gene_id + '_' + gene_name
            amp_seq = eles[11]

            exon_starts = eles[8].split(',')
            exon_ends = eles[9].split(',')
            exon_count = len(exon_starts) - 1
            seq_overlap_list = []
            coding_seq = ''

            for i in range(0, exon_count):
                exon_id = genome_gene + '_exon_' + exon_starts[i] + '_' + exon_ends[i]
                exon_seq = exonid2fa[exon_id]
                seq_overlap = overlap_seq_finder(amp_seq, exon_seq)
                seq_overlap_list.append(seq_overlap)

            coding_seq = max(seq_overlap_list, key = len)
            newline = line + '\t' + coding_seq + '\n'
            newlines.append(newline)

    with open(os.path.join(dir_path, master_file2out), 'w') as master2out:
        master2out.write(header + '\tcoding_seq\n')
        for newline in newlines:
            master2out.write(newline)

#################################################################################################################################

    
