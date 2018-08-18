
import numpy as np
import pandas as pd
import glob, yaml, json, xmltodict
import IPython

files_dict = json.load(open("metadata/gdc_metadata.json"))

file_uuid_map = {}
for entry in files_dict:
	file_uuid_map[entry["file_id"]] = entry

files_per_case = {}
case_project = {}

for dir_name in glob.glob("files/*"):
	uuid = dir_name[6:]
	entry = file_uuid_map[uuid]
	case_id = entry ['cases'][0]['case_id']
	project_id = entry['cases'][0]["project"]["project_id"]
	data_file = dir_name + "/" + entry ['file_name']
	data_type = entry['data_type']

	if case_id not in files_per_case:
		files_per_case[case_id] = [(data_file, data_type)]
	else:
		files_per_case[case_id].append((data_file, data_type))
	case_project[case_id] = project_id

	print (uuid, files_per_case[case_id])
	print ("\n\n")


data = {}
for case in files_per_case.keys():

	old_id = ""
	case_data = {}
	for data_file, data_type in files_per_case[case]:
		data_file = str(data_file)
		if data_type == 'Clinical Supplement':
			case_data['clinical_data_file'] = data_file
			ind = data_file.find("TCGA")
			old_id = data_file[ind:ind+12]
		if data_type == 'Gene Expression Quantification':
			case_data['gene_expression_file'] = data_file
		if data_type == 'miRNA Expression Quantification':
			case_data['mirna_expression_file'] = data_file
		if data_type == 'Isoform Expression Quantification':
			case_data['isoform_expression_file'] = data_file

	files = ["clinical_data_file", 'gene_expression_file', "mirna_expression_file", 'isoform_expression_file']
	if "clinical_data_file" in case_data:
		print (case_data)
		case_data['project'] = str(case_project[case])
		data[old_id] = case_data

yaml.dump(data, open("processed/case_files_locs.yaml",'w'))
