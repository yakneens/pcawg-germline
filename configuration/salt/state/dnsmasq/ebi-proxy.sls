/etc/dnsmasq.d/10-ebi-proxy:    
  file.append:
    - text: "server=192.168.0.5"
    - makedirs: True
