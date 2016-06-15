from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from sqlalchemy import or_, and_
import os.path
import datetime
import uuid
import json


def write_config_to_file(config, config_location):
    
    run_uuid = str(uuid.uuid4())
    
    my_file = open("{}/{}.json".format(config_location, run_uuid), "w")
    json.dump(config, my_file)
    my_file.close()

results_path = "/shared/data/results/regenotype_pcawg_tumors/"
tissue_type = "tumor"
config_location = "/tmp/missing_configs"

counter = 0

contig_names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                "20", "21", "22", "X", "Y"]

Base = automap_base()        
sample_engine = create_engine('postgresql://pcawg_admin:pcawg@postgresql.service.consul:5432/pcawg_sample_tracking')        
Base.prepare(sample_engine, reflect=True)        
            
PCAWGSample = Base.classes.pcawg_samples
SampleLocation = Base.classes.sample_locations
 
sample_session = Session(sample_engine)

if tissue_type == "normal":
    sample_id = PCAWGSample.normal_wgs_alignment_gnos_id
    sample_location = SampleLocation.normal_sample_location
else:
    sample_id = PCAWGSample.tumor_wgs_alignment_gnos_id
    sample_location = SampleLocation.tumor_sample_location 

flag = True
    
for root, dirs, files in os.walk(results_path):
    
    if flag:
        flag = False
        continue
    
    missing_contigs = []
    this_sample_id = os.path.basename(root)
    
    for my_contig in contig_names:
        if not os.path.isfile(root + "/" + this_sample_id + "_" + my_contig + ".vcf.gz"):
            missing_contigs.append(my_contig)
    

    if len(missing_contigs) > 0:

        current_sample = sample_session.query(PCAWGSample.index.label("index"), sample_id.label("sample_id"), sample_location.label("sample_location")).\
            join(SampleLocation, PCAWGSample.index == SampleLocation.donor_index).\
            filter(sample_id == this_sample_id).first()
        
        this_config_data = {"sample": {
                                    "donor_index": current_sample.index,
                                    "sample_id": str.split(current_sample.sample_id, ",")[0],
                                    "sample_location": current_sample.sample_location,
                                },
                            "contig_whitelist": missing_contigs
                            }
            
        write_config_to_file(this_config_data, config_location)
        

