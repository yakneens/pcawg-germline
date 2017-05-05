shared_filesystems:
  file.append:
    - name: /etc/dnsmasq.d/10-germline
    - text: "server=/em-isi-3104.ebi.ac.uk/10.35.104.201"
    - makedirs: True