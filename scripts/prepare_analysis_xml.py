import xml.etree.ElementTree as ET

import uuid
import os.path
import datetime
import re
from xml.etree.ElementTree import SubElement
from subprocess import check_output

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


def prepare_sample_files_for_submission(submission_base_path, sample_id, sample_files):
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
        new_file_full_path = "{}/{}/{}".format(submission_base_path, new_submission_id)
        copy_command = "cp {}/{} {}".format(sample_path_prefix,
                                            sample_file,
                                            new_filename)
        
        call_command(copy_command, "copy")
        
        md5 = check_output("md5sum " + new_file_full_path)
        new_files.append((new_filename, md5))
        
    return new_files, new_submission_id

def prepare_submission_metadata_from_template(metadata_template_location, sample_id, submission_full_path, dcc_project_code, submitter_donor_id, icgc_or_tcga):
    analysis_template = ET.parse(metadata_template_location)

    analysis_template.find("./ANALYSIS")[0].set("analysis_date", str(datetime.datetime.now()))
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
    
prepared_files, new_submission_id = prepare_sample_files_for_submission(submission_base_path, sample_id, sample_files)    
prepare_submission_metadata_from_template(metadata_template_location, sample_id, submission_base_path + "/" + new_submission_id, dcc_project_code, submitter_donor_id, icgc_or_tcga)



