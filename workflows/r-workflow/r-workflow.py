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


def run_R(**kwargs):
    "results_base_path": "/shared/data/results/r_density_bins_by_project/",
    "results_local_path": "/tmp/r_density_bins_by_project/",
    "donor_meta_path": "/shared/data/density_bins_input/donor_meta.RData",
    "deletion_ranges_path": "/shared/data/density_bins_input/deletion_ranges.RData",
    "snv_ranges_path": "/shared/data/density_bins_input/snv_ranges.RData",
    "deletion_carrier_mask_path": "shared/data/density_bins_input/deletion_carrier_mask.RData"

    config = get_config(kwargs)
    
    
    
    result_path = config["results_local_path"]
    results_base_path = config["results_base_path"]
    donor_meta_path = config["donor_meta_path"]
    deletion_ranges_path = config["deletion_ranges_path"]
    snv_ranges_path = config["snv_ranges_path"]
    deletion_carrier_mask_path = config["deletion_carrier_mask_path"]
    
    
    
    if (not os.path.isdir(result_path)):
        logger.info(
            "Results directory {} not present, creating.".format(result_path))
        os.makedirs(result_path_prefix)

    
    bcftools_command = 'bcftools {} {} {} -o {}'.\
        format(bcftools_operation,
               bcftools_flags,
               filenames_string,
               result_filename)
    

    call_command(bcftools_command, "bcftools")

    generate_tabix(result_filename, config)
    copy_result(result_filename, sample_id, config)
   
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

dag = DAG("R", default_args=default_args,
          schedule_interval=None, concurrency=20000, max_active_runs=20000)


start_analysis_run_task = PythonOperator(
    task_id="start_analysis_run",
    python_callable=start_analysis_run,
    provide_context=True,
    dag=dag)

r_task = PythonOperator(
    task_id="R",
    python_callable=run_R,
    provide_context=True,
    dag=dag)

r_task.set_upstream(start_analysis_run_task)

complete_analysis_run_task = PythonOperator(
    task_id="complete_analysis_run",
    python_callable=complete_analysis_run,
    provide_context=True,
    dag=dag)

complete_analysis_run_task.set_upstream(r_task)
