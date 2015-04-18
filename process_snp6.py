s = open("/Users/siakhnin/Downloads/SNP6_uuid_barcode.csv", "r")
w = open("/Users/siakhnin/Downloads/WGS_donor_uuid_barcode.csv", "r")
v = open("/Users/siakhnin/Downloads/PCAWG_donor_uuid_vcf_uuid.csv", "r")
f = open("/Users/siakhnin/Downloads/snp6_vcf_tcga_uuids.txt","r")

w_input = w.readlines()
w_input = w_input[0].split("\r")
w_input.pop(0)
wgs_donor_ids = {}
for e in w_input:
   el = e.split(",")
   wgs_donor_ids[el[1]] = el[0]
    
vcf_file_ids = {}    
v_input = v.readlines()
v_input = v_input[0].split("\r")
for e in v_input:
	el = e.split(",")
	vcf_file_ids[el[0]] = el[1]

s_input = s.readlines()

s_input = s_input[0].split("\r")
s_input.pop(0)
result_list = []


final_snp6_uuids = {}
f_input = f.readlines()
for e in f_input:
	final_snp6_uuids[e.strip()] = 1

o = open("/Users/siakhnin/Downloads/SNP6_to_vcf_file_uuid_map.csv","w")
i =0
for e in s_input:
    cur = e.split(",")
    barcode = cur[1]
    donor_barcode = barcode[0:12]
    if final_snp6_uuids.has_key(cur[0]) and wgs_donor_ids.has_key(donor_barcode) and vcf_file_ids.has_key(wgs_donor_ids[donor_barcode]):
    	result_list.append(cur[0] + "," + vcf_file_ids[wgs_donor_ids[donor_barcode]] + "," + cur[4])
    	i = i + 1
    	#o.write(cur[0] + "," + vcf_file_ids[wgs_donor_ids[donor_barcode]] + "\n")
for i in result_list:
	print i