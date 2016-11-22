#!/usr/bin/env python
import os
import uuid
import json
import argparse


def parse_args():
    my_parser = argparse.ArgumentParser()
    
    sub_parsers = my_parser.add_subparsers()
    
    create_configs_parser = sub_parsers.add_parser("create-configs", conflict_handler='resolve')
    create_configs_parser.add_argument("-c", "--config_location", help="Path to a directory where the generated config files will be stored.", dest="config_location", required=True)
    create_configs_parser.set_defaults(func=create_configs_command)
    
    my_args = my_parser.parse_args()
    
    return my_args

def write_config_to_file(config, config_location):
    
    run_uuid = str(uuid.uuid4())
    
    my_file = open("{}/{}.json".format(config_location, run_uuid), "w")
    json.dump(config, my_file)
    my_file.close()

def generate_config_objects(chromosome_names, config_location):
    for chromosome_name in chromosome_names:
    
        
        this_config_data = {"r_param": {
                                "-o": chromosome_name,
                                }
                            }
        
        yield this_config_data
        

def create_configs_command(args):

    config_location = args.config_location
    
    chromosome_names = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    
    if (not os.path.isdir(config_location)):
        os.makedirs(config_location)
        
    for config in generate_config_objects(chromosome_names, config_location):
        write_config_to_file(config, config_location)
    
if __name__ == '__main__':
    args = parse_args()
    args.func(args)
    

    
