#!/usr/bin/env python
import os
import uuid
import json
import argparse

from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from sqlalchemy import or_, and_
import re



def parse_args():
    my_parser = argparse.ArgumentParser()

    sub_parsers = my_parser.add_subparsers()

    create_configs_parser = sub_parsers.add_parser(
        "create-configs", conflict_handler='resolve')
    create_configs_parser.add_argument(
        "-n", "--num_runs", help="Number of runs to create configurations for.", dest="num_runs", required=False, type=int)
    create_configs_parser.add_argument(
        "-c", "--config_location", help="Path to a directory where the generated config files will be stored.", dest="config_location", required=True)
    create_configs_parser.add_argument(
        "-s", "--sample_location", help="Path to a directory where the sample files are stored.", dest="sample_location", required=True)
    
    create_configs_parser.set_defaults(func=create_configs_command)

    my_args = my_parser.parse_args()

    return my_args


def write_config_to_file(config, config_location):

    run_uuid = str(uuid.uuid4())

    my_file = open("{}/{}.json".format(config_location, run_uuid), "w")
    json.dump(config, my_file)
    my_file.close()


def create_configs_command(args):

    num_runs = args.num_runs
    config_location = args.config_location
    sample_location = args.sample_location
    
    num_configs = 0
    
    Base = automap_base()        
    sample_engine = create_engine('postgresql://pcawg_admin:pcawg@postgresql.service.consul:5432/pcawg_sample_tracking')        
    Base.prepare(sample_engine, reflect=True)        
                
    PCAWGSample = Base.classes.pcawg_samples
                
    sample_session = Session(sample_engine)
    
    sample_id = PCAWGSample.normal_wgs_aliquot_id
    
    
            
    
    if (not os.path.isdir(config_location)):
        os.makedirs(config_location)

    
    for root, dirs, files in os.walk(sample_location):
        if num_runs == None or num_configs < num_runs:
            sample_uuid = os.path.basename(root)
            sample_data = sample_session.query(PCAWGSample.dcc_project_code.label("dcc_project_code"), PCAWGSample.submitter_donor_id.label("submitter_donor_id")).\
                            filter(sample_id == sample_uuid).\
                            first()
            icgc_or_tcga = "ICGC"
            if (re.search("_US$", sample_data.dcc_project_code)):
                icgc_or_tcga = "TCGA"
                               
                            
            this_config_data = {
                                "sample": {
                                           "sample_id": sample_uuid,
                                           "path_prefix": root,
                                           "sample_files": files
                                           },
                                "dcc_project_code": sample_data.dcc_project_code,
                                "submitter_donor_id": sample_data.submitter_donor_id,
                                "icgc_or_tcga": icgc_or_tcga
                                }
            write_config_to_file(this_config_data, config_location)
            num_configs = num_configs + 1
        elif num_configs >= num_runs:
            return


    sample_session.close()
    sample_engine.dispose()
    

if __name__ == '__main__':
    args = parse_args()
    args.func(args)
