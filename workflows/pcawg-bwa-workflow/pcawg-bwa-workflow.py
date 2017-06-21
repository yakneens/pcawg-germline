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

from subprocess import check_output


def run_bwa(**kwargs):
    config = get_config(kwargs)
    sample = get_sample(kwargs)
    
    sample_id = sample["sample_id"]
    sample_path_prefix = sample["path_prefix"]
    sample_files = sample["sample_files"]
    dcc_project_code = config["dcc_project_code"]
    submitter_donor_id = config["submitter_donor_id"]
    
    
    
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
