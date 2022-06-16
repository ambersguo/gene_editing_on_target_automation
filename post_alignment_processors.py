##################################################################################################################################
#                                                     post_alignment_processors
# report_maker
# 
# 
#
#
#
##################################################################################################################################

import os, glob, re, string
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
from collections import defaultdict
from functions import guide_master_parser

########################################################report_maker#############################################################

def report_maker (scripts_path, mode, run_filter, gid2batchline, gid2pg_setting, gid2ew_setting):

    path_batchfile = scripts_path.replace('scripts', 'fastq')
    path_input_q = scripts_path.replace('scripts', 'fastq/CRISPRessoBatch_on_batchfile')
    path_input_s = scripts_path.replace('scripts', 'fastq_s30/CRISPRessoBatch_on_batchfile')
    path_input_e = scripts_path.replace('scripts', 'fastq_ew_s30/CRISPRessoBatch_on_batchfile')
    path_input_p = scripts_path.replace('scripts', 'fastq_pf/CRISPRessoBatch_on_batchfile')
    path2out = scripts_path.replace('scripts', 'processing')
    path_inputs = []
    if 'q' in run_filter:
        path_inputs.append(path_input_q)
    if 's' in run_filter:
        path_inputs.append(path_input_s)
    if 'e' in run_filter:
        path_inputs.append(path_input_e)
    if 'p' in run_filter:
        path_inputs.append(path_input_p)

    # parse sample names
    sample_list_all = []
    output2sample = {}
    with open(os.path.join(path_batchfile, 'batchfile.tsv'), 'r') as batch2in:
        header = batch2in.readline()
        for line in batch2in:
            line = line.rstrip('\n')
            eles = re.split(r'\t',line)
            sample = eles[0]
            sample_list_all.append(sample)
            sample_output = sample.replace('-', '')
            output2sample[sample_output] = sample

    gid2es = {}
    ampid2es = {}
    for gid, ew_setting in gid2ew_setting.items():
        gid2es[gid] = ew_setting[-1]

    for gid, batchline in gid2batchline.items():
        batch_eles = batchline.split('\t')
        amp_seq = batch_eles[2]
        amp_id = batch_eles[3]
        es = gid2es[gid]
        ampid2es[amp_id] = es

    # make general final reports
    for path_input in path_inputs:
        if 'fastq/CRISPRessoBatch' in path_input:
            group_input = 'q30'
        elif 'fastq_s30/CRISPRessoBatch' in path_input:
            group_input = 's30'
        elif 'fastq_ew_s30/CRISPRessoBatch' in path_input:
            group_input = 'ews30'
        elif 'fastq_pf/CRISPRessoBatch' in path_input:
            group_input = 'pf'

        samples_output = []
        indel_ratios = []
        sub_ratios = []
        newlines_report = []
        sample_output2fs = {}

        with open(os.path.join(path_input, 'CRISPRessoBatch_quantification_of_frameshift_splicing.txt'), 'r') as file2in_fs:
            header = file2in_fs.readline()
            for line in file2in_fs:
                line = line.rstrip('\n')
                eles = line.split('\t')
                sample_output = eles[0]
                reads_fs = int(eles[3])
                sample_output2fs[sample_output] = reads_fs

        with open(os.path.join(path_input, 'CRISPRessoBatch_quantification_of_editing_frequency.txt'), 'r') as file2in:
            header = file2in.readline()
            for line in file2in:
                line = line.rstrip('\n')
                eles = line.split('\t')
                sample_output = eles[0]
                sample = output2sample[sample_output]
                samples_output.append(sample)
                sample_eles = sample.split('-')
                
                if 'G0004' in eles[3]:
                    gid =  'GEx-G0004'
                    plate_well_id = 'pl5'
                else:
                    gid = 'GEx-' + eles[3] + 'CC'
                    plate_well_id = eles[6]
                
                batchline = gid2batchline[gid]
                batch_eles = batchline.split('\t')
                amp_id = batch_eles[3]

                reads_input = int(eles[4])
                reads_aligned = int(eles[6])
                reads_ins = int(eles[10])
                reads_del = int(eles[11])
                reads_indel = reads_ins + reads_del - int(eles[16])
                reads_sub = int(eles[12])
                reads_fs = sample_output2fs[sample_output]
                aligned_ratio = '{:.2f}'.format((float(reads_aligned)/reads_input)*100)
                ins_ratio = '{:.2f}'.format((float(reads_ins)/reads_aligned)*100)
                del_ratio = '{:.2f}'.format((float(reads_del)/reads_aligned)*100)
                indel_ratio = '{:.2f}'.format((float(reads_indel)/reads_aligned)*100)
                sub_ratio = '{:.2f}'.format((float(eles[12])/reads_aligned)*100)
                fs_ratio = '{:.2f}'.format((float(reads_fs)/reads_aligned)*100)
                indel_ratios.append(float(indel_ratio))
                sub_ratios.append(float(sub_ratio))

                QC_mark = ''
                if reads_input < 1000:
                    QC_mark += 'low_raw_read_count_'
                if float(aligned_ratio) < 10:
                    QC_mark += 'low_percentage_alignment_'
                if float(sub_ratio) > 5:
                    QC_mark += 'high_percentage_substitution_'

                newline_report = sample + '\t' + gid + '\t' + amp_id + '\t' + \
                                 str(reads_input) + '\t' + str(reads_aligned) + '\t' + str(reads_ins)  + '\t' + str(reads_del)  + '\t' + str(reads_indel) + '\t' + str(reads_sub) + '\t' + \
                                 str(aligned_ratio) + '\t' + str(ins_ratio) + '\t' + str(del_ratio) + '\t' + str(indel_ratio) + '\t' + str(sub_ratio) + '\t' + str(QC_mark)

                newlines_report.append(newline_report)

        for sample in sample_list_all:
            if sample not in samples_output:
                newline_report = sample + '\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tpre_alignment_QC_failed'
                newlines_report.append(newline_report)

        df = pd.DataFrame([newline.split('\t') for newline in newlines_report])
        df[0] = df[0].astype(str)
        df_sorted = df.sort_values(by=[0], ascending = [True])
        newlines_report_sorted = df_sorted.values.tolist()

        with open(os.path.join(path2out, group_input + '_final_report.txt'), 'w') as file2report:
            file2report.write('Sample\tGuide_id\tAmplicon_id\tReads_input\tReads_aligned\tReads_insertion\tReads_deletion\tReads_indel\tReads_substitution' +  \
                              '\tAligned_percentage(%)\tInsertion_percentage(%)\tDeletion_percentage(%)\tIndel_percentage(%)\tSubstitution_percentage(%)\tQC_filter\n')

            for newline_report_sorted in newlines_report_sorted:
                newline = '\t'.join([str(ele) for ele in newline_report_sorted])
                file2report.write(newline + '\n')

        print('Completed: General final report generation for ' + mode + ': ' + group_input)

    # make base-editing specific reports
        if 'BE' in mode:
            for filename in os.listdir(path_input):
                if 'Nucleotide_percentage_summary.txt' in filename:
                    filename_eles = filename.split('.')
                    amp_id = filename_eles[0]
                    es = ampid2es[amp_id]

                    subs = []
                    sample_sub2ratio = defaultdict(float)
                    sample2es_sum = defaultdict(float)
                    amp_samples = []
                    with open(os.path.join(path_input, filename), 'r') as befile2in:
                        header = befile2in.readline()
                        header = header.rstrip('\n')
                        header_eles = header.split('\t')
                        es_col = int(es) + 1
                        amp_es = header_eles[es_col]

                        for line in befile2in:
                            line = line.rstrip('\n')
                            eles = line.split('\t')
                            sample_output = eles[0]
                            sample = output2sample[sample_output]
                            amp_samples.append(sample)
                            sample_es = eles[1]
                            sample_es_ratio = float(eles[es_col]) * 100

                            if sample_es != 'N' and sample_es != '-' and sample_es != amp_es:
                                sub = amp_es + '-to-' + sample_es
                                subs.append(sub)
                                sample_sub = sample + '&' + sub
                                sample_sub2ratio[sample_sub] = sample_es_ratio
                                sample2es_sum[sample] += sample_es_ratio
                
                    amp_samples_uniq = list(set(amp_samples))
                    amp_samples_uniq.sort()
                    subs_uniq = list(set(subs))
                    subs_uniq.sort()

                    befile2out_name = group_input + '_' + amp_id + "_base_editing_frequency_report.txt"
                    with open(os.path.join(path2out, befile2out_name), 'w') as befile2out:
                        header2out = 'Sample\tAmplicon_id\tSum_base_editing_percentage_at_editing_site(%)\t' + subs_uniq[0] + '(%)\t' + subs_uniq[1] + '(%)\t' + subs_uniq[2] + '(%)\n'
                        befile2out.write(header2out)
                        for amp_sample in amp_samples_uniq:
                            sample_es_ratio_1 = sample_es_ratio_2 = sample_es_ratio_3 = 0
                            for sub in subs_uniq:
                                sample_sub = amp_sample + '&' + sub
                                sample_es_ratio = sample_sub2ratio[sample_sub]
                                if sub == subs_uniq[0]:
                                    sample_es_ratio_1 = float(sample_es_ratio)
                                elif sub == subs_uniq[1]:
                                    sample_es_ratio_2 = float(sample_es_ratio)
                                elif sub == subs_uniq[2]:
                                    sample_es_ratio_3 = float(sample_es_ratio)

                            newline = amp_sample + '\t' + amp_id + '\t' + str(sample2es_sum[amp_sample]) + '\t' + str(sample_es_ratio_1) + '\t' + str(sample_es_ratio_2) + '\t' + str(sample_es_ratio_3)
                            befile2out.write(newline + '\n')
        print('Base editing mode frequency report generation completed for:' + group_input)



##################################################################################################################################


