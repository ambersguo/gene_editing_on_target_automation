############################################################################################################################################
#                                                     functions
# setting_parser
# run_info_parser
# guide_master_parser
# fastq_parser
# sample_name_parser
# check_machine
# fastq_idhash
# seq_fast_generator
# overlap_seq_finder
#
############################################################################################################################################



import os, sys, gzip, re, string
from difflib import SequenceMatcher
from collections import defaultdict

##################################################################setting_parser#############################################################

def setting_parser (setting_filename):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, setting_filename), 'r') as setting_file:
        setting_all = {}
        for line in setting_file:
            line = line.rstrip('\n')
            eles = re.split(r'\t',line)
            setting_all[eles[0]] = eles[1] + '&' + eles[2]
    
    setting4batch_line = ""
    fprimerid2seq = {}
    rprimerid2seq = {}
    fguideid2seq = {}
    rguideid2seq = {}
    w_leftF = w_rightF = w_leftR = w_rightR = offset = 0

    for key in setting_all:
        if 'run' in key:
            run_id, rest = setting_all[key].split('&')
        elif 'primer' in key:
            primer_id, primer_seq = setting_all[key].split('&')
            if 'f' in setting_all[key]:
                fprimerid2seq[primer_id] = primer_seq
            elif 'r' in setting_all[key]:
                rprimerid2seq[primer_id] = primer_seq
        elif 'guide' in key:
            guide_id, guide_seq = setting_all[key].split('&')
            setting4batch_line = guide_seq + '\t' + guide_id + '\t'
            fguideid2seq[guide_id] = guide_seq
            guide_r = guide_seq[::-1]
            guide_temp = guide_r.translate(str.maketrans("AT", "TA"))
            guide_rc = guide_temp.translate(str.maketrans("CG", "GC"))
            rguideid2seq[guide_id] = guide_rc
        elif 'ref' in key:
            ref_id, ref_seq = setting_all[key].split('&')
            setting4batch_line += ref_seq + '\t' + ref_id + '\t'
            if guide_seq in ref_seq:
                guide_mark = 'guide_forward'
            else:
                guide_mark = 'guide_reverse'
        elif 'offset' in key:
            off_id, offset = setting_all[key].split('&')           
        elif 'editing_window' in key:
            ewid, ew = setting_all[key].split('&')
            setting4batch_line +=  '30\t' + ew + '\t1'
            w_leftF, w_rightF = ew.split('-')
            w_leftF = int(w_leftF) + int(offset)
            w_rightF = int(w_rightF) + int(offset)
            w_leftR = len(ref_seq) - w_rightF + int(offset)*2
            w_rightR = len(ref_seq) - w_leftF + int(offset)*2

    setting_pg = [fprimerid2seq, rprimerid2seq, fguideid2seq, rguideid2seq]
    setting_es = [w_leftF, w_rightF, w_leftR, w_rightR]
    settings = [run_id, setting4batch_line, setting_pg, setting_es, guide_mark]
    return settings

#############################################################################################################################################





##################################################################run_info_parser############################################################

def run_info_parser (run_info_filename):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, run_info_filename), 'r') as run_info_file:
        run_info_file.readline()
        id2info = {}
        for line in run_info_file:
            line = line.rstrip('\n')
            eles = line.split('\t')
            id2info[eles[0]] = eles[1]
    return id2info

#############################################################################################################################################





##################################################################guide_master_parser########################################################

def guide_master_parser (guide_master_filename, window_range, min_paired):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, guide_master_filename), 'r') as master2in:
        master2in.readline()

        gid2batchline = {}
        gid2pg_setting = {}
        gid2ew_setting = {}
        edit_site = 0
        settings = []
        for line in master2in:
            line = line.rstrip('\n')
            eles = line.split('\t')
            gid = eles[0]
            guide_name = eles[1]
            guide_seq = eles[2]
            amp_id = eles[3]
            amp_seq = eles[4]
            primerF = eles[5]
            primerR = eles[6]
            mode = eles[7]

            guide_r = guide_seq[::-1]
            guide_temp = guide_r.translate(str.maketrans("AT", "TA"))
            guide_rc = guide_temp.translate(str.maketrans("CG", "GC"))

            if guide_seq in amp_seq:
                amp_seq_pre, amp_seq_after = amp_seq.split(guide_seq)
                edit_site = len(amp_seq_pre) + len(guide_seq) - 3
                guide_mark = 'guide_forward'
            else:
                if guide_rc in amp_seq:
                    amp_seq_pre, amp_seq_after = amp_seq.split(guide_rc)
                    edit_site = len(amp_seq_pre) + 3
                    guide_mark = 'guide_reverse'
                else:
                    sys.exit("Error: guide sequence not found in amplicon sequence, please check\n")
            
            edit_window_left = edit_site - int(window_range)
            edit_window_right = edit_site + int(window_range)
            if edit_window_left < 0:
                sys.exit("Error: editing window range is too large, please check\n")
            if edit_window_right > len(amp_seq):
                sys.exit("Error: editing window range is too large, please check\n")
            edit_window = str(edit_window_left) + '-' + str(edit_window_right)

            gid2batchline[gid] =  guide_seq + '\t' + gid + '\t' + amp_seq + '\t' + amp_id + '\t' + '30\t' + edit_window + '\t' + min_paired

            gid2pg_setting[gid] = [primerF, primerR, guide_seq, guide_rc, guide_mark]

            w_leftF = edit_window_left - 10
            w_rightF = edit_window_left + 10
            w_leftR = len(amp_seq) - w_rightF - 10
            w_rightR = len(amp_seq) - w_leftF + 10
            gid2ew_setting[gid] = [w_leftF, w_rightF, w_leftR, w_rightR, edit_site]

        settings = [gid2batchline, gid2pg_setting, gid2ew_setting]
        print("Guide master data parsing completed.")

        return settings

#############################################################################################################################################






##################################################################fastq_parser###############################################################

def fastq_parser (fastq_path, sample_divider):
    samples = []
    fastq_parsed = []
    sample2R1 = {}
    sample2R2 = {}
    R1_count = R2_count = 0
    for filename in os.listdir(fastq_path):
        if '.fastq.gz' in filename:
            sample, rest = filename.split(sample_divider)
            samples.append(sample)
            if '_R1' in filename:
                sample2R1[sample] = filename
                R1_count += 1
            elif '_R2' in filename:
                sample2R2[sample] = filename
                R2_count += 1
    fastq_parsed = [samples, sample2R1, R1_count, sample2R2, R2_count]
    return fastq_parsed

#############################################################################################################################################





########################################################sample_name_parser##################################################################

def sample_name_parser(sample_name):
    sample_name_info = []
    eles = sample_name.split('-')
    seq_id = eles[1]
    plate_well_id = eles[4]
    gid =  'GEx-' + eles[3]
    sample_name_info = [seq_id, plate_well_id, gid]
    return sample_name_info

#############################################################################################################################################






########################################################check_machine########################################################################

def check_machine (fastq_path):
    for filename in os.listdir(fastq_path):
        if '.fastq.gz' in filename:
            with gzip.open(os.path.join(fastq_path, filename), 'rt') as file2in:
                for line in file2in:
                    line = str(line)
                    machine = line.split(':')[0]
                    return machine

#############################################################################################################################################






########################################################fastq_idhash#########################################################################

def fastq_idhash (fastq_path, machine):
    idhash_list = []
    Fid2seq = {}
    Rid2seq = {}
    Fid2qual = {}
    Rid2qual = {}
    mark = ''
    for filename in os.listdir(fastq_path):
        if 'fastq.gz' in filename:
            if '_R1_' in filename:
                mark = 'R1'
            elif '_R2_' in filename:
                mark = 'R2'

            with gzip.open(os.path.join(fastq_path, filename), 'rt') as file2in:
                line_count = 1
                for line in file2in:
                    line = str(line)
                    line = line.rstrip('\n')
                    if line.startswith(machine) and line_count == 1:
                        seqid, rest = line.split(' ')
                        line_count += 1
                    elif line_count == 2:
                        seq = line
                        if mark == 'R1':
                            Fid2seq[seqid] = seq
                        elif mark == 'R2':
                            Rid2seq[seqid] = seq
                        line_count += 1
                    elif line_count == 3:
                        line_count += 1
                    elif line_count == 4:
                        qual = line
                        if mark == 'R1':
                            Fid2qual[seqid] = qual
                        elif mark == 'R2':
                            Rid2qual[seqid] = qual
                        line_count = 1

    idhash_list = [Fid2seq, Rid2seq, Fid2qual, Rid2qual]
    return idhash_list

#############################################################################################################################################






########################################################seq_fast_generator########################################################

def seq_fast_generator(refpath, chrom, seqstart, seqend):

    with open(refpath + chrom + '.fa', 'r') as ref:
        if seqend < seqstart:
            raise Exception('Incorrect coordinates!\n')
        header = ref.readline()
        header = header.rstrip('\n')
        if not re.match(r'^>', header):
            raise Exception('The reference file is not in fasta format!\n')
        
        headeroffset = len(header)
        lineref = ref.readline()
        basesperline = len(lineref) - 1

        lfstostart = int((seqstart/basesperline))
        lfswithin = int(seqend/basesperline) - lfstostart
        startbyte = headeroffset + lfstostart + seqstart
        bytestoread = seqend + lfswithin - seqstart
        
        ref.seek(startbyte, 0)
        seq = ref.read(bytestoread)
        seq = re.sub('\s', '', seq)
        return seq

#################################################################################################################################







########################################################overlap_seq_finder#######################################################

def overlap_seq_finder(S1, S2):
  M = [[0]*(1+len(S2)) for i in range(1+len(S1))]
  longest, x_longest = 0, 0
  for x in range(1,1+len(S1)):
    for y in range(1,1+len(S2)):
        if S1[x-1] == S2[y-1]:
            M[x][y] = M[x-1][y-1] + 1
            if M[x][y]>longest:
                longest = M[x][y]
                x_longest  = x
        else:
            M[x][y] = 0
  return S1[x_longest-longest: x_longest]

#################################################################################################################################


