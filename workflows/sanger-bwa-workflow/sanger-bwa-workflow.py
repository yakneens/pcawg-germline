from airflow import DAG
from airflow.operators import BashOperator, PythonOperator
from datetime import datetime, timedelta

import os
import logging
from subprocess import call

from string import join

import tracker.model
from tracker.model.analysis_run import *
from tracker.util.workflow_common import *

import uuid
import os.path
import datetime
import re
import json

from subprocess import check_output

def get_cwl_config_template():
    fp = open("config-template.json")
    cwl_config = json.load(fp)
    fp.close()
    return cwl_config

def prepare_cwl_config(cwl_config, config, sample):
    reference_location = config["reference_location"]
    bwa_index_location = config["bwa_index_location"]
    generate_cram = config["generate_cram"]
    output_mapping = config["output_mapping"]
    #output_mapping = {"out_bam":".bam", "out_bai":".bam.bai", "out_bas":".bam.bas", "out_md5":".bam.md5", "out_met":".bam.met", "out_maptime":".bam.maptime"}
    
    sample_id = sample["sample_id"]
    sample_location = sample["sample_location"]
    
    result_path_prefix = "/tmp/"
    
    if (not os.path.isdir(result_path_prefix)):
        logger.info(
            "Results directory {} not present, creating.".format(result_path_prefix))
        os.makedirs(result_path_prefix)
    
    cwl_config["reference"]["path"] = reference_location
    cwl_config["bwa_idx"]["path"] = bwa_index_location
    
    cwl_config["bams_in"][0]["path"] = sample_location
    
    cwl_config["cram"] = generate_cram
    
    for out_key in output_mapping:
        cwl_config[out_key]["path"] = "{}{}{}{}".format(result_path_prefix, sample_id, "_mapped", output_mapping[out_key])  
    
    
    cwl_config_location = config["cwl_config_location"]
    
    fp = open(cwl_config_location,"w")
    json.dump(cwl_config, fp)
    fp.close()
    
    return cwl_config_location
    

def run_bwa(**kwargs):
    config = get_config(kwargs)
    sample = get_sample(kwargs)
    
    sample_id = sample["sample_id"]
    sample_location = sample["sample_location"]
    reference_location = config["reference_location"]
    sample_path_prefix = sample["path_prefix"]
    dcc_project_code = config["dcc_project_code"]
    submitter_donor_id = config["submitter_donor_id"]
    
    cwl_flags = config["cwl_flags"]
    
    result_path_prefix = config["results_local_path"] + "/" + sample_id
    
    cwl_config_location = prepare_cwl_config(get_cwl_config_template(), config, sample)
    
    
    cwl_file_location = config["cwl_file_location"]
    
    cwl_command = "{} {} --outdir {} {} {}".\
        format("cwl-runner",
               cwl_flags,
               result_path_prefix,
               cwl_file_location,
               cwl_config_location)
    call_command(cwl_command, "cwl")

    #copy_result(compressed_sample_filename, sample_id, config)
    #"cwl-runner --debug --non-strict --enable-dev Dockstore.cwl examples/sample_configs.local.json"
    
    return submission_sample_location



   
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'start_date': datetime.datetime(2020, 01, 01),
    'email': ['airflow@airflow.com'],
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}

dag = DAG("pcawg_bwa", default_args=default_args,
          schedule_interval=None, concurrency=50, max_active_runs=50)


start_analysis_run_task = PythonOperator(
    task_id="start_analysis_run",
    python_callable=start_analysis_run,
    provide_context=True,
    dag=dag)

run_bwa_task = PythonOperator(
    task_id="run_bwa",
    python_callable=run_bwa,
    provide_context=True,
    dag=dag)

run_bwa_task.set_upstream(start_analysis_run_task)

complete_analysis_run_task = PythonOperator(
    task_id="complete_analysis_run",
    python_callable=complete_analysis_run,
    provide_context=True,
    dag=dag)

complete_analysis_run_task.set_upstream(run_bwa_task)
