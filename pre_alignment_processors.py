##################################################################################################################################
#                                                     pre_alignment_processor
# bathfile_maker
# primer_filter
# editing_window_s30_filter
#
#
#
##################################################################################################################################




import os, gzip, re, string
import pandas as pd
from collections import defaultdict
from functions import setting_parser, guide_master_parser, fastq_parser, sample_name_parser, check_machine, fastq_idhash




########################################################batchfile_maker###########################################################

def batchfile_maker (fastq_path, sample_divider, gid2batchline, mode):

    fastq_parsed = fastq_parser(fastq_path, sample_divider)
    [samples, sample2R1, R1_count, sample2R2, R2_count] = fastq_parsed
    newlines_batchfile = []
    samples = list(set(samples))
    samples.sort()
    if R1_count:
        if R1_count == R2_count:
            for sample in samples:
                gid = sample_name_parser(sample)[-1]
                filename_R1 = sample2R1[sample]
                filename_R2 = sample2R2[sample]
                batch_line = sample + '\t' + filename_R1 + '\t' + filename_R2 + '\t' + gid2batchline[gid]
                newlines_batchfile.append(batch_line)
                print('fastq files found as pair-end reads: ' + filename_R1 + ', ' + filename_R2 + '\n')
        elif R2_count == 0:
            for sample in samples:
                gid = sample_name_parser(sample)[-1]
                filename_R1 = sample2R1[sample]
                batch_line = sample + '\t' + filename_R1 + '\t' + gid2batchline[gid]
                newlines_batchfile.append(batch_line)
                print('fastq files found as single-end reads: ' + filename_R1 + '\n')
    else:
        print('Error: no fastq.gz file found in your directory, please check\n')

    with open(os.path.join(fastq_path, 'batchfile.tsv'), 'w') as file2out:
        if R1_count == R2_count:
            header = 'name\tfastq_r1\tfastq_r2\tguide_seq\tguide_name\tamplicon_seq\tamplicon_name\tq\tqwc\tmin_paired_end_reads_overlap\n'
        elif R2_count == 0:
            header = 'name\tfastq_r1\tguide_seq\tguide_name\tamplicon_seq\tamplicon_name\tq\tqwc\tmin_paired_end_reads_overlap\n'
        file2out.write(header)
        for newline in newlines_batchfile:
            file2out.write(newline + '\n')
    
    message = 'Completed: making batchfile for mode ' + mode + ' in folder ' + fastq_path + '\n'
    print(message)

#############################################################################################################################################





#######################################################primer_filter########################################################################

def primer_filter (fastq_path, fastq_pf_path, stats_path, gid2pg_setting, fastq_hash, machine):

    print('Started: Primer filtering.\n')

    filename2out = 'stats_primer_seq_counts_pre_alignment.txt'
    newlines = []
    [Fid2seq, Rid2seq, Fid2qual, Rid2qual] = fastq_hash
    
    for filename in os.listdir(fastq_path):
        primer_seqids = []
        primerF = ''
        primerR = ''
        primer_inread = ''
        mark = ''
        if 'fastq.gz' in filename:
            filename_eles = filename.split('_')
            sample_eles = filename_eles[0].split('-')
            gid = sample_eles[2] + '-' + sample_eles[3]
            primerF = gid2pg_setting[gid][0]
            primerR = gid2pg_setting[gid][1]

            if '_R1_' in filename:
                sample = filename_eles[0] + '_pf_R1'
                primer_inread = primerF
                mark = 'R1'
            elif '_R2_' in filename:
                sample = filename_eles[0] + '_pf_R2'
                primer_inread = primerR
                mark = 'R2'

            with gzip.open(os.path.join(fastq_path, filename), 'rt') as file2in:
                check = total_count = primer_count = 0

                for line in file2in:
                    line = str(line)
                    line = line.rstrip('\n')
  
                    if line.startswith(machine):
                        seqid, rest = line.split(' ')
                        check = 1
                        total_count += 1
                    elif check:
                        seq = line
                        check = 0

                        if mark == 'R1':
                            if primerF in seq:
                                primer_seqids.append(seqid)
                                primer_count += 1
                        elif mark == 'R2':
                            if primerR in seq:
                                primer_seqids.append(seqid)
                                primer_count += 1

                newline = sample + '\t' + str(total_count) + '\t' + str(primer_count) + '\t' + str('{:.2f}'.format((float(primer_count/total_count)*100)))
                newlines.append(newline)

            fastq2out_name = sample + '.fastq'
            with open(os.path.join(fastq_pf_path, fastq2out_name), 'w') as fastq2out:
                if mark == 'R1':
                    for seqid in primer_seqids:
                        seq = Fid2seq[seqid]
                        qual = Fid2qual[seqid]
                        fastq_line = seqid + '\n' + seq + '\n+\n' + qual + '\n'
                        fastq2out.write(fastq_line)
                elif mark == 'R2':
                    for seqid in primer_seqids:
                        seq = Rid2seq[seqid]
                        qual = Rid2qual[seqid]                     
                        fastq_line = seqid + '\n' + seq + '\n+\n' + qual + '\n'
                        fastq2out.write(fastq_line)

        df = pd.DataFrame([newline.split('\t') for newline in newlines])
        df[0] = df[0].astype(str)
        df_sorted = df.sort_values(by=[0], ascending = [True])
        newlines_sorted = df_sorted.values.tolist()

    with open(os.path.join(stats_path, filename2out), 'w') as file2out:
        header = 'sample_id\ttotal_read_count\tread_w_primer_count\tread_w_primer_percentage(%)\n'
        file2out.write(header)
        for newline_sorted in newlines_sorted:
            newline = '\t'.join([str(ele) for ele in newline_sorted])
            file2out.write(newline + '\n')
    print('Completed: Primer filtering.\n')

#############################################################################################################################################





########################################################editing_window_s30_filter###########################################################

def ew_filter (fastq_path, fastq_ew_path, gid2ew_setting, fastq_hash, score):

    print('Started: Editing window filtering.')

    [Fid2seq, Rid2seq, Fid2qual, Rid2qual] = fastq_hash

    samples = []
    for filename in os.listdir(fastq_path):
        if 'fastq.gz' in filename:
            filename_eles = filename.split('_')
            sample = filename_eles[0]
            samples.append(sample)

    samples = list(set(samples))
    for sample in samples:
        seqids = []
        seqidFs = []
        seqidRs = []
        filenames = []

        gid = sample_name_parser(sample)[-1]
        [w_leftF, w_rightF, w_leftR, w_rightR, es] = gid2ew_setting[gid]

        for filename in os.listdir(fastq_path):
            sample2 = sample + '_'
            if sample2 in filename:
                filenames.append(filename)
        filenames.sort()

        for filename in filenames:
            if 'fastq.gz' in filename:
                if '_R1' in filename:
                    with gzip.open(os.path.join(fastq_path, filename), 'rt') as file2inF:
                        line_count = 1
                        for line in file2inF:
                            line = str(line)
                            line = line.rstrip('\n')
                            if line_count == 1:
                                seqid = line.split(' ')[0]
                            elif line_count == 4:
                                qualf = line
                                line_count = 0
                                qf_pass = 1
                                wf_qual = qualf[w_leftF : w_rightF+1]
                                wf_qual_list = list(wf_qual)
                                for q_asc in wf_qual_list:
                                    q = ord(q_asc) - 33
                                    if q < int(score):
                                        qf_pass = 0
                                if qf_pass:
                                    seqidFs.append(seqid)
                            line_count += 1

                elif '_R2' in filename:
                    with gzip.open(os.path.join(fastq_path, filename), 'rt') as file2inR:
                        line_count = 1
                        for line in file2inR:
                            line = str(line)
                            line = line.rstrip('\n')
                            if line_count == 1:
                                seqid = line.split(' ')[0]
                            elif line_count == 4:
                                qualr = line
                                line_count = 0
                                qr_pass = 1
                                if w_rightR - len(qualr) > -1:
                                    wr_qual = qualr[w_leftR :]
                                else:
                                    wr_qual = qualr[w_leftR : w_rightR+1]
                                wr_qual_list = list(wr_qual)
                                for q_asc in wr_qual_list:
                                    q = ord(q_asc) - 33
                                    if q < int(score):
                                        qr_pass = 0
                                if qr_pass:
                                    if seqid in seqidFs:
                                        seqidRs.append(seqid)
                            line_count += 1

        filename2outFfq = sample + '_ew_s30_R1.fastq'
        filename2outRfq = sample + '_ew_s30_R2.fastq'

        with open(os.path.join(fastq_ew_path, filename2outFfq), 'w') as file2outFfq:
            for seqid in seqidRs:
                newlineF_fq = seqid + '\n' + Fid2seq[seqid] + '\n' + '+\n' + Fid2qual[seqid] + '\n'
                file2outFfq.write(newlineF_fq)

        with open(os.path.join(fastq_ew_path, filename2outRfq), 'w') as file2outRfq:
            for seqid in seqidRs:
                newlineR_fq = seqid + '\n' + Rid2seq[seqid] + '\n' + '+\n' + Rid2qual[seqid] + '\n'
                file2outRfq.write(newlineR_fq)

    print('Completed: Editing window filtering.')

#############################################################################################################################################



