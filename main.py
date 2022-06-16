################################################################################################################
#                                gene_editing on-target evaluation pipeline
#
#
#
#
################################################################################################################

import os, sys, gzip, re, string
import subprocess as sp
from contextlib import contextmanager
from functions import run_info_parser, guide_master_parser, fastq_parser, check_machine, fastq_idhash
from pre_alignment_processors import batchfile_maker, primer_filter, ew_filter
from post_alignment_processors import report_maker


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def main():

    # run info parsing
    id2info = run_info_parser('run_info.txt')
    run_id = id2info['run_id']
    mode = id2info['mode']
    ew_score = id2info['ew_score']
    run_filter = id2info['run_filter']
    window_range = id2info['window_range']
    min_paired = id2info['min_paired']

    # guide master parsing
    [gid2batchline, gid2pg_setting, gid2ew_setting] = guide_master_parser('guide_master.txt', window_range, min_paired)

    # info of working directories
    run_mkdir_cmd = "mkdir " + run_id
    run_path = '../' + run_id
    run_scripts_path = run_path + '/scripts'
    run_fastq_path = run_path + '/fastq'
    run_fastq_s30_path = run_path + '/fastq_s30'
    run_fastq_pf_path = run_path + '/fastq_pf'
    run_fastq_ew_path = run_path + '/fastq_ew_s30'
    run_fastq_qc_path = run_path + '/fastq_qc'
    run_processing_path = run_path + '/processing'

    # info of html links
    run_aws_folder_link_q30 = 's3://s3-informatics/shared-data/experiments/gene_editing_ngs/ambersguo/' + run_id + '/'
    run_aws_folder_link_s30 = run_aws_folder_link_q30 + 's30/'
    run_aws_folder_link_pf = run_aws_folder_link_q30 + 'pf/'
    run_aws_folder_link_ew = run_aws_folder_link_q30 + 'ews30/'
    run_link_head = 'https://s3.modernatx.net/s3-informatics/shared-data/experiments/gene_editing_ngs/ambersguo/' 
    run_aws_folder_sync_q30_cmd = 'aws s3 sync CRISPRessoBatch_on_batchfile ' + run_aws_folder_link_q30 + 'CRISPRessoBatch_on_batchfile'
    run_aws_html_sync_q30_cmd = 'aws s3 cp CRISPRessoBatch_on_batchfile.html ' + run_aws_folder_link_q30
    run_aws_html_q30_link = run_link_head + run_id + '/CRISPRessoBatch_on_batchfile.html'
    run_aws_folder_sync_s30_cmd = 'aws s3 sync CRISPRessoBatch_on_batchfile ' + run_aws_folder_link_s30 + 'CRISPRessoBatch_on_batchfile'
    run_aws_html_sync_s30_cmd = 'aws s3 cp CRISPRessoBatch_on_batchfile.html ' + run_aws_folder_link_s30
    run_aws_html_s30_link = run_link_head + run_id + '/s30/CRISPRessoBatch_on_batchfile.html'
    run_aws_folder_sync_pf_cmd = 'aws s3 sync CRISPRessoBatch_on_batchfile ' + run_aws_folder_link_pf + 'CRISPRessoBatch_on_batchfile'
    run_aws_html_sync_pf_cmd = 'aws s3 cp CRISPRessoBatch_on_batchfile.html ' + run_aws_folder_link_pf
    run_aws_html_pf_link = run_link_head + run_id + '/pf/CRISPRessoBatch_on_batchfile.html'
    run_aws_folder_sync_ew_cmd = 'aws s3 sync CRISPRessoBatch_on_batchfile ' + run_aws_folder_link_ew + 'CRISPRessoBatch_on_batchfile'
    run_aws_html_sync_ew_cmd = 'aws s3 cp CRISPRessoBatch_on_batchfile.html ' + run_aws_folder_link_ew
    run_aws_html_ew_link = run_link_head + run_id + '/ews30/CRISPRessoBatch_on_batchfile.html'
    run_aws_folder_sync_qc_cmd = 'aws s3 sync multiqc_data ' + run_aws_folder_link_q30 + 'multiqc'
    run_aws_html_sync_qc_cmd = 'aws s3 cp multiqc_report.html ' + run_aws_folder_link_q30
    run_aws_html_qc_link = run_link_head + run_id + '/multiqc_report.html'

    docker_fastqc_cmd = "docker run --rm \
    -u 2048:2048 \
    -v ${PWD}:/data \
    --workdir /data \
    465747649825.dkr.ecr.us-east-1.amazonaws.com/compsci/fastqc:1.0.0 \
    fastqc *.fastq.gz"

    docker_multiqc_cmd = "docker run --rm \
    -u 2048:2048 \
    -v ${PWD}:/data \
    --workdir /data \
    --interactive \
    465747649825.dkr.ecr.us-east-1.amazonaws.com/compsci/multiqc:358bd98339170d908653b450ace8e4afd623a706 \
    multiqc ."

    docker_crispresso_q30_cmd = "docker run --rm \
    -u 2048:2048 \
    -v ${PWD}:/data \
    --workdir /data \
    465747649825.dkr.ecr.us-east-1.amazonaws.com/compsci/crispresso2:1.1.0 \
    CRISPRessoBatch \
    --batch_settings batchfile.tsv \
    --use_legacy_insertion_quantification \
    --n_processes 12 \
    --skip_failed"

    docker_crispresso_s30_cmd = docker_crispresso_q30_cmd + " -s 30"

    docker_crispresso_BE_q30_cmd = docker_crispresso_q30_cmd + " --base_editor_output"

    docker_crispresso_BE_s30_cmd = docker_crispresso_BE_q30_cmd + " -s 30"

    with cd('../'):
        sp.call(run_mkdir_cmd, shell = True)

    with cd(run_path):
        sp.call("mkdir fastq fastq_s30 fastq_pf fastq_ew_s30 fastq_qc processing scripts", shell = True)
        sp.call("cp ../scripts_new_ew_ns3/* ./scripts", shell = True)

    with cd(run_fastq_qc_path):
        sp.call("aws s3 cp s3://mtx-molbio/NGS_Scott/" + run_id + "/ ./ --recursive", shell = True)
        sp.call(docker_fastqc_cmd, shell = True)
        sp.call(docker_multiqc_cmd, shell = True)

    with cd(run_fastq_path):
        sp.call("aws s3 cp s3://mtx-molbio/NGS_Scott/" + run_id + "/ ./ --recursive", shell = True)

    # generate batchfile for q30 filter
    with cd(run_scripts_path):
        run_dir_path = os.path.dirname(os.path.realpath(__file__))
        full_fastq_path = run_dir_path.replace('scripts', 'fastq')
        batchfile_maker(full_fastq_path, '_S', gid2batchline, mode)

    # run with q30 filter
    if 'q' in run_filter:
        with cd(run_fastq_path):
            print('Started: on-target evaluation with q30 filter for mode ' + mode + '.\n')
            if mode == 'DSB':
                sp.call(docker_crispresso_q30_cmd, shell = True)
            elif mode == 'BE':
                sp.call(docker_crispresso_BE_q30_cmd, shell = True)
            print('Completed: on-target evaluation with q30 filter for mode ' + mode + '.\n')

    # run with s30 filter
    if 's' in run_filter:
        with cd(run_fastq_s30_path):
            print('Started: on-target evaluation with s30 filter for mode ' + mode + '.\n')
            sp.call("aws s3 cp s3://mtx-molbio/NGS_Scott/" + run_id + "/ ./ --recursive", shell = True)
            sp.call("cp ../fastq/batchfile.tsv ./", shell = True)
            if mode == 'DSB':
                sp.call(docker_crispresso_s30_cmd, shell = True)
            elif mode == 'BE':
                sp.call(docker_crispresso_BE_s30_cmd, shell = True)
            print('Completed: on-target evaluation with s30 filter for mode ' + mode + '.\n')

    # make fastq hash
    print('Started: making fastq hash.\n')
    machine = check_machine(full_fastq_path)
    fastq_hash = fastq_idhash(full_fastq_path, machine)
    print('Completed: making fastq hash.\n')

    # run with ews30 filter
    if 'e' in run_filter:
        with cd(run_fastq_ew_path):
            print('Started: on-target evaluation with editing window filter for mode ' + mode + '.\n')
            full_fastq_ew_path = os.path.dirname(os.path.realpath(__file__))
            full_fastq_path = full_fastq_path.replace('fastq_ew_s30', 'fastq')
            ew_filter(full_fastq_path, full_fastq_ew_path, gid2ew_setting, fastq_hash, ew_score)
            sp.call("gzip *.fastq", shell = True)
            batchfile_maker(full_fastq_ew_path, '_ew', gid2batchline, mode)
            if mode == 'DSB':
                sp.call(docker_crispresso_q30_cmd, shell = True)
            elif mode == 'BE':
                sp.call(docker_crispresso_BE_q30_cmd, shell = True)
            print('Completed: on-target evaluation with editing window filter for mode ' + mode + '.\n')

    # run with primer filter
    if 'p' in run_filter:
        with cd(run_fastq_pf_path):
            print('Started: on-target evaluation with primer filter for mode ' + mode + '.\n')
            full_fastq_pf_path = os.path.dirname(os.path.realpath(__file__))
            full_fastq_path = full_fastq_pf_path.replace('fastq_pf', 'fastq')
            full_stats_path = full_fastq_pf_path.replace('fastq_pf', 'processing')
            primer_filter(full_fastq_path, full_fastq_pf_path, full_stats_path, gid2pg_setting, fastq_hash, machine)
            sp.call("gzip *.fastq", shell = True)
            batchfile_maker(full_fastq_pf_path, '_pf', gid2batchline, mode)
            sp.call(docker_crispresso_q30_cmd, shell = True)
            sp.call(run_aws_folder_sync_pf_cmd, shell = True)
            sp.call(run_aws_html_sync_pf_cmd, shell = True)
            print('Completed: on-target evaluation with primer filter for mode ' + mode + '.\n')
    
    # generate reports
    with cd(run_scripts_path):
        run_dir_path = os.path.dirname(os.path.realpath(__file__))
        report_maker(run_dir_path, mode, run_filter, gid2batchline, gid2pg_setting, gid2ew_setting)

    # push html outputs to aws s3
    with cd(run_fastq_qc_path):
        sp.call(run_aws_folder_sync_qc_cmd, shell = True)
        sp.call(run_aws_html_sync_qc_cmd, shell = True)

    if 'q' in run_filter:
        with cd(run_fastq_path):
            sp.call(run_aws_folder_sync_q30_cmd, shell = True)
            sp.call(run_aws_html_sync_q30_cmd, shell = True)

    if 's' in run_filter:
        with cd(run_fastq_s30_path):
            sp.call(run_aws_folder_sync_s30_cmd, shell = True)
            sp.call(run_aws_html_sync_s30_cmd, shell = True)

    if 'e' in run_filter:
        with cd(run_fastq_ew_path):
            sp.call(run_aws_folder_sync_ew_cmd, shell = True)
            sp.call(run_aws_html_sync_ew_cmd, shell = True)
    
    if 'p' in run_filter:
        with cd(run_fastq_pf_path):
            sp.call(run_aws_folder_sync_pf_cmd, shell = True)
            sp.call(run_aws_html_sync_pf_cmd, shell = True)

    with cd(run_processing_path):
        full_run_processing_path = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(full_run_processing_path, 'output_htmls.txt'), 'w') as file2out:
            file2out.write(run_aws_html_qc_link + '\n')
            if 'q' in run_filter:
                file2out.write(run_aws_html_q30_link + '\n')
            if 's' in run_filter:
                file2out.write(run_aws_html_s30_link + '\n')
            if 'e' in run_filter:
                file2out.write(run_aws_html_ew_link + '\n')
            if 'p' in run_filter:
                file2out.write(run_aws_html_pf_link + '\n')

            print('Completed: pushing output htmls to aws s3 bucket.\n')



if __name__ == '__main__':
    main()
