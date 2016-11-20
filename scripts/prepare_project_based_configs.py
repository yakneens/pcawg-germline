#!/usr/bin/env python
import sys
import os
import uuid
from time import sleep
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from sqlalchemy import or_, and_
from sqlalchemy import distinct
import json
import logging
import argparse

from tracker.util import connection

def parse_args():
    my_parser = argparse.ArgumentParser()
    
    sub_parsers = my_parser.add_subparsers()
    
    create_configs_parser = sub_parsers.add_parser("create-configs", conflict_handler='resolve')
    create_configs_parser.add_argument("-c", "--config_location", help="Path to a directory where the generated config files will be stored.", dest="config_location", required=True)
    create_configs_parser.set_defaults(func=create_configs_command)
    
    my_args = my_parser.parse_args()
    
    return my_args

def get_project_names():
    
    #PCAWG Samples are in their own database
    Base = automap_base()        
    sample_engine = create_engine('postgresql://pcawg_admin:pcawg@postgresql.service.consul:5432/pcawg_sample_tracking')        
    Base.prepare(sample_engine, reflect=True)        
                
    PCAWGSample = Base.classes.pcawg_samples
                
    sample_session = Session(sample_engine)
    
    project_names = sample_session.query(distinct(PCAWGSample.dcc_project_code)).all()
        
    sample_session.close()
    sample_engine.dispose()
    
    
    
    return project_names

def write_config_to_file(config, config_location):
    
    run_uuid = str(uuid.uuid4())
    
    my_file = open("{}/{}.json".format(config_location, run_uuid), "w")
    json.dump(config, my_file)
    my_file.close()

def generate_config_objects(project_names, config_location):
    for project_name in project_names:
    
        
        this_config_data = {"r_param": {
                                "-p": project_name,
                                }
                            }
        
        yield this_config_data
        

def create_configs_command(args):

    config_location = args.config_location
    
    project_names = get_project_names()
    
    if (not os.path.isdir(config_location)):
        os.makedirs(config_location)
        
    for config in generate_config_objects(project_names, config_location):
        write_config_to_file(config, config_location)
    
if __name__ == '__main__':
    args = parse_args()
    args.func(args)
    

    
