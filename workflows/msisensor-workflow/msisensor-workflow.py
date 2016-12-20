from airflow import DAG
from airflow.operators import BashOperator, PythonOperator
from datetime import datetime, timedelta

import os
import logging
from subprocess import call

import tracker.model
#from tracker.model.analysis_run import *
#from tracker.util.workflow_common import *

def run_msisensor(**kwargs):
    config = get_config(kwargs)
    logger.debug("Config - {}".format(config))

    sample = get_sample(kwargs)
    sample_id = sample["tumor_wgs_aliquot_id"]
    tumor_location = sample["tumor_bam_path"]
    normal_location = sample["normal_bam_path"]

    result_path_prefix = config["results_local_path"] + "/" + sample_id
    if (not os.path.isdir(result_path_prefix)):
        logger.info(
            "Results directory {} not present, creating.".format(result_path_prefix))
        os.makedirs(result_path_prefix)

    result_filename_base = "{}/{}".format(result_path_prefix, sample_id)
    # MSIsensor creates four output files:
    # base, base_somatic, base_germline, base_dis

    msisensor_path = config["msisensor"]["path"]
    msisensor_flags = config["msisensor"]["flags"]
    microsatellites_ref = config["microsatellites_ref"]
    n_threads = 1

    msisensor_cmd = "{path} -d {ref} -n {normal} -t {tumor} -o {output} -b {threads} {flags}".\
        format(path = msisensor_path, ref = microsatellites_ref,
               normal = normal_location, tumor = tumor_location, output = result_filename_base,
               threads = n_threads, flags = msisensor_flags)

    call_command(msisensor_cmd, "msisensor")
    call_command('rm {}_dis'.format(result_filename_base))
    results = []
    results.append(compress_sample(result_filename_base, config))
    results.append(compress_sample(result_filename_base + '_germline', config))
    results.append(result_filename_base + '_somatic')
    for res_fn in results:
        copy_result(res_fn, sample_id, config)

default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'start_date': datetime(2020, 01, 01),
    'email': ['airflow@airflow.com'],
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}

dag = DAG("msisensor", default_args = default_args,
          schedule_interval = None, concurrency = 10000,
          max_active_runs = 2000)

start_analysis_run_task = PythonOperator(
    task_id="start_analysis_run",
    python_callable=start_analysis_run,
    provide_context=True,
    dag=dag)

msisensor_task = PythonOperator(
    task_id = 'msisensor',
    python_callable = run_msisensor,
    provide_context = True,
    dag = dag)

msisensor_task.set_upstream(start_analysis_run_task)

complete_analysis_run_task = PythonOperator(
    task_id = "complete_analysis_run",
    python_callable = complete_analysis_run_task,
    provide_context = True,
    dag = dag)

complete_analysis_run_task.set_upstream(msisensor_task)
