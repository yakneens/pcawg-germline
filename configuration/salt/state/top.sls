base:
  '*':
    - consul
    - dnsmasq
    - dnsmasq.ebi-proxy
    - collectd
    - elastic.filebeat
    - elastic.packetbeat
    - ntp
  'G@roles:consul-bootstrap':
    - consul.bootstrap
  'G@roles:consul-server':
    - consul.server
  'G@roles:consul-client':
    - consul.client
  'G@roles:monitoring-server':
    - influxdb
    - grafana 
  'G@roles:worker':
    - dnsmasq.germline-share
    - celery
    - airflow
    - airflow.load-workflows
    - airflow.worker
    - butler.tracker
    - cwltool
    - docker 
  'G@roles:tracker':
    - dnsmasq.germline-share
    - run-tracking-db.set_db_url
    - celery
    - airflow
    - airflow.init-db
    - airflow.load-workflows
    - airflow.server
    - jsonmerge
    - butler.tracker    
  'G@roles:db-server':
    - postgres
    - run-tracking-db
    - run-tracking-db.create_tables
    - grafana.createdb
    - airflow.airflow-db
    - sample-tracking-db
  'G@roles:elasticsearch':
    - elastic.search
    - elastic.logstash
    - elastic.kibana
    - celery
  'G@roles:tmp_germline':
    - biotools.freebayes
    - biotools.htslib
    - biotools.samtools
    - biotools.delly
  'G@roles:job-queue':
    - rabbitmq

