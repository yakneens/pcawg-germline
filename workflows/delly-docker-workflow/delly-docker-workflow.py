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
import distutils.util

from subprocess import check_output

def get_cwl_config_template():
    fp = open("/opt/airflow/dags/pcawg-germline/delly-docker-workflow/config-template.json")
    cwl_config = json.load(fp)
    fp.close()
    return cwl_config

def prepare_cwl_config(cwl_config, config, sample):
    reference_location = config["reference_location"]
    exclude_regions_location = config["exclude_regions_location"]
    reference_gc_location = config["reference_gc_location"]
    gc_wig_location = config["gc_wig_location"]
    sv_pon_location = config["sv_pon_location"]
    num_cores = config["num_cores"]
    
    normal_sample_id = sample["normal_sample_id"]
    normal_sample_location = sample["normal_sample_location"]
    
    tumour_sample_id = sample["tumour_sample_id"]
    tumour_sample_location = sample["tumour_sample_location"]
    
    result_path_prefix = "/tmp/"
    
    if (not os.path.isdir(result_path_prefix)):
        logger.info(
            "Results directory {} not present, creating.".format(result_path_prefix))
        os.makedirs(result_path_prefix)
    
    cwl_config["reference-gz"]["path"] = reference_location
    cwl_config["exclude-reg"]["path"] = exclude_regions_location
    cwl_config["reference-gc"]["path"] = reference_gc_location
    cwl_config["gc_wig"]["path"] = gc_wig_location
    cwl_config["sv-collection"]["path"] = sv_pon_location
    cwl_config["normal-bam"]["path"] = normal_sample_location
    cwl_config["tumor-bam"]["path"] = tumour_sample_location
    
    cwl_config["ncpu"] = num_cores

    
    fp = open(cwl_config_location,"w")
    json.dump(cwl_config, fp)
    fp.close()
    
    return cwl_config_location
    

def run_delly(**kwargs):
    config = get_config(kwargs)
    sample = get_sample(kwargs)
    
    normal_sample_id = sample["normal_sample_id"]
    normal_sample_location = sample["normal_sample_location"]
    
    tumour_sample_id = sample["tumour_sample_id"]
    tumour_sample_location = sample["tumour_sample_location"]
    
    cwl_flags = config["cwl_flags"]
    
    result_path_prefix = config["results_base_path"] + "/" + normal_sample_id
    
    cwl_config_location = prepare_cwl_config(get_cwl_config_template(), config, sample)
    
    
    cwl_file_location = config["cwl_file_location"]
    
    cwl_command = "{} {} --outdir {} --tmp-outdir-prefix /pan-prostate/docker_temp/ {} {}".\
        format("cwl-runner",
               cwl_flags,
               result_path_prefix,
               cwl_file_location,
               cwl_config_location)
    call_command("id", "id")
    call_command(cwl_command, "cwl")

   
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

dag = DAG("docker_delly", default_args=default_args,
          schedule_interval=None, concurrency=500, max_active_runs=500)


start_analysis_run_task = PythonOperator(
    task_id="start_analysis_run",
    python_callable=start_analysis_run,
    provide_context=True,
    dag=dag)

run_delly_task = PythonOperator(
    task_id="run_delly",
    python_callable=run_delly,
    provide_context=True,
    dag=dag)

run_delly_task.set_upstream(start_analysis_run_task)

complete_analysis_run_task = PythonOperator(
    task_id="complete_analysis_run",
    python_callable=complete_analysis_run,
    provide_context=True,
    dag=dag)

complete_analysis_run_task.set_upstream(run_bwa_task)
