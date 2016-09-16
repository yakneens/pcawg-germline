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

import xml.etree.ElementTree as ET

import uuid
import os.path
import datetime
import re
from xml.etree.ElementTree import SubElement
from subprocess import check_output

def prepare_sample_files_for_submission(sample_path_prefix, submission_base_path, sample_id, sample_files):
    if (not os.path.isdir(submission_base_path)):
        logger.info(
            "Submission directory {} not present, creating.".format(submission_base_path))
        os.makedirs(submission_base_path)
            
    new_submission_id = str(uuid.uuid4())
    os.makedirs("{}/{}".format(submission_base_path, new_submission_id))
    
    new_files = []
    
    for sample_file in sample_files:
        ext = sample_file[sample_file.find("."):]
        new_filename = "{}.butler-freebayes-1-0-0.{}.germline.snv{}".format(sample_id, datetime.datetime.now().strftime('%Y%m%d'), ext)
        new_file_full_path = "{}/{}/{}".format(submission_base_path, new_submission_id, new_filename)
        copy_command = "cp {}/{} {}".format(sample_path_prefix,
                                            sample_file,
                                            new_file_full_path)
        
        call_command(copy_command, "copy")
        logger.info(new_file_full_path)
        md5 = check_output(["md5sum", new_file_full_path]).split(" ")[0]
        new_files.append((new_filename, md5))
        
    return new_files, new_submission_id

def prepare_submission_metadata_from_template(prepared_files, metadata_template_location, 
                                              sample_id, submission_full_path, dcc_project_code, 
                                              submitter_donor_id, icgc_or_tcga):
    analysis_template = ET.parse(metadata_template_location)

    analysis_template.find("./ANALYSIS").set("analysis_date", datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
    seq_labels = analysis_template.findall("./ANALYSIS/ANALYSIS_TYPE/REFERENCE_ALIGNMENT/SEQ_LABELS/SEQUENCE")
    
    for label in seq_labels:
        label.set("data_block_name", sample_id)
        
    analysis_template.find("./ANALYSIS/TARGETS/TARGET").set("refname", sample_id)
    
    analysis_template.find("./ANALYSIS/DATA_BLOCK").set("name", sample_id)
    
    files_el = analysis_template.find("./ANALYSIS/DATA_BLOCK/FILES")
    
    for new_file in prepared_files:
        file_el = SubElement(files_el, "FILE")
        file_el.set("filename", new_file[0])
        
        el_type = ""
        
        if(re.search(".vcf.gz.tbi", new_file[0]) != None):
            el_type = "idx"
        else:
            el_type = "vcf"
        
        file_el.set("filetype", el_type)
        file_el.set("checksum_method", "MD5")
        file_el.set("checksum", new_file[1])
    
    
    
    analysis_template.find("./ANALYSIS/ANALYSIS_ATTRIBUTES/ANALYSIS_ATTRIBUTE[TAG='dcc_project_code']/VALUE").text = dcc_project_code
    analysis_template.find("./ANALYSIS/ANALYSIS_ATTRIBUTES/ANALYSIS_ATTRIBUTE[TAG='submitter_donor_id']/VALUE").text = submitter_donor_id
    
    analysis_template.write(submission_full_path + "/analysis.xml", "utf-8")      
    

def prepare_metadata(**kwargs):
    config = get_config(kwargs)
    sample = get_sample(kwargs)
    
    sample_id = sample["sample_id"]
    sample_path_prefix = sample["path_prefix"]
    sample_files = sample["sample_files"]
    dcc_project_code = config["dcc_project_code"]
    submitter_donor_id = config["submitter_donor_id"]
    icgc_or_tcga = config["icgc_or_tcga"]
    
    metadata_template_location = config["metadata_template_location"]
    submission_base_path = config["submission_base_path"]
    
    prepared_files, new_submission_id = prepare_sample_files_for_submission(sample_path_prefix, 
                                                                            submission_base_path, 
                                                                            sample_id, 
                                                                            sample_files)
    
    submission_sample_location =  submission_base_path + "/" + new_submission_id
    
    prepare_submission_metadata_from_template(prepared_files, metadata_template_location, 
                                              sample_id, submission_sample_location, 
                                              dcc_project_code, submitter_donor_id, icgc_or_tcga)
    
    return submission_sample_location

def submit_metadata(**kwargs):
    submission_sample_location = kwargs["ti"].xcom_pull(task_ids="prepare_metadata")
    
    config = get_config(kwargs)
    sample = get_sample(kwargs)
    icgc_or_tcga = config["icgc_or_tcga"]
    
    gnos = config["gnos"]
    destination_repo_mapping = config["destination_repo_mapping"]
    
    selected_repo = gnos[destination_repo_mapping[icgc_or_tcga]]
    
    cgsubmit_command = "cgsubmit -u {} -s {} -c {}".format(submission_sample_location, 
                                                           selected_repo["url"], 
                                                           selected_repo["key_location"])
    call_command(cgsubmit_command, "cgsubmit")

def upload_sample(**kwargs):
    submission_sample_location = kwargs["ti"].xcom_pull(task_ids="prepare_metadata")
    
    config = get_config(kwargs)
    sample = get_sample(kwargs)
    icgc_or_tcga = config["icgc_or_tcga"]
    
    gnos = config["gnos"]
    destination_repo_mapping = config["destination_repo_mapping"]
    
    selected_repo = gnos[destination_repo_mapping[icgc_or_tcga]]
    
    gtupload_command = "/opt/gtupload/cghub/bin/gtupload {} -c {} -v".format(submission_sample_location + "/manifest.xml", 
                                                     selected_repo["key_location"])
    call_command(gtupload_command, "gtupload", cwd=submission_sample_location)
    

   
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

dag = DAG("gtupload", default_args=default_args,
          schedule_interval=None, concurrency=20000, max_active_runs=50)


start_analysis_run_task = PythonOperator(
    task_id="start_analysis_run",
    python_callable=start_analysis_run,
    provide_context=True,
    dag=dag)

metadata_task = PythonOperator(
    task_id="prepare_metadata",
    python_callable=prepare_metadata,
    provide_context=True,
    dag=dag)

metadata_task.set_upstream(start_analysis_run_task)

cgsubmit_task = PythonOperator(
    task_id="submit_metadata",
    python_callable=submit_metadata,
    provide_context=True,
    dag=dag)

cgsubmit_task.set_upstream(metadata_task)

gtupload_task = PythonOperator(
    task_id="upload_sample",
    python_callable=upload_sample,
    provide_context=True,
    dag=dag)

gtupload_task.set_upstream(cgsubmit_task)

complete_analysis_run_task = PythonOperator(
    task_id="complete_analysis_run",
    python_callable=complete_analysis_run,
    provide_context=True,
    dag=dag)

complete_analysis_run_task.set_upstream(gtupload_task)
