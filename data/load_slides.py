
import numpy as np
import pandas as pd
import os, glob, yaml, json, shutil, xmltodict

import IPython

data = yaml.load(open("processed/case_files_locs.yaml"))
df = pd.read_csv("manifests/gdc_tissue_slides_manifest.tsv", sep='\t')
prefix = "cloud-data/NCI-GDC/legacy/TCGA/"
cache_dir = "/media/nik/Seagate Backup Plus Drive/tissue-slides/"


for case in sorted(data.keys()):
	project = data[case]['project']
	directory = prefix + project + "/Other/Tissue_slide_image/"

	case_files_df = df[df['filename'].str.contains(case)]
	ids = case_files_df['id'].tolist()
	filenames = case_files_df['filename'].tolist()

	data[case]['slides'] = []
	for i, (file_id, filename) in enumerate(zip(ids, filenames)):
		data_file = directory + file_id + "/" + filename
		
		try:
			cache_file = cache_dir + case + "_slide" + str(i) + ".svs"
			print ("Caching to: ", cache_file)
			shutil.copy(data_file, cache_file)
			data_file = cache_file
		except:
			pass

		data[case]['slides'].append(data_file)

		print (data_file)

yaml.dump(data, open("processed/case_files_locs.yaml",'w'))

