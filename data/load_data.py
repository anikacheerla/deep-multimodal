
import numpy as np
import pandas as pd
import glob, yaml, json, xmltodict

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.cross_validation import cross_val_predict, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
from sklearn.metrics import roc_auc_score, accuracy_score

import IPython

def load_data_from_path(xml_dict, *args):
	for spec in args:
		try:
			key_matching = [key for key in xml_dict.keys() if spec in key][0]
		except: return False

		xml_dict = xml_dict[key_matching]
	return xml_dict


data = yaml.load(open("processed/case_files_locs.yaml"))
print len(data)
clinical_data = []
clinical_data2 = []

race_lookup = []
drug_lookup = []
disease_lookup = []

keys = ["recurrence", "vital_status", "gender", "race", "age", "disease", "mirna_data", "gene_data"] 


for case in data.keys():
	entry = data[case]
	print (entry['clinical_data_file'])
	dd = xmltodict.parse(open(entry['clinical_data_file']))

	clin = load_data_from_path(dd, "tcga", "patient", "new_tumor_events", "new_tumor_event_after_initial_treatment", "text")
	if clin: data[case]['recurrence'] = 1 if clin == "YES" else 0

	clin = load_data_from_path(dd, "tcga", "patient", 'vital_status', 'text')
	if clin: data[case]['vital_status'] = 1 if clin == "Alive" else 0

	clin = load_data_from_path(dd, "tcga", "patient", 'gender', 'text')
	if clin: data[case]['gender'] = 1 if clin == "MALE" else 0

	clin = load_data_from_path(dd, "tcga", "patient", 'race_list', 'race', 'text')
	clin2 = load_data_from_path(dd, "tcga", "patient", 'ethnicity', 'text')
	if clin and clin2:
		if clin2 == "HISPANIC OR LATINO":
			clin = "LATINO"
		if "NATIVE" in clin:
			clin = "ASIAN"

		if clin not in race_lookup: race_lookup.append(clin)
		data[case]['race'] = race_lookup.index(clin)

	clin = load_data_from_path(dd, "tcga", "patient", 'age_at_initial', 'text')
	if clin: data[case]['age'] = int(clin)

	clin = load_data_from_path(dd, "tcga", "patient", 'drugs', 'drug', "drug_name", 'text')
	if clin:
		if clin not in drug_lookup: drug_lookup.append(clin)
		data[case]['drug'] = drug_lookup.index(clin)

	clin = load_data_from_path(dd, "tcga", "patient", 'histologic_grade', 'text')
	
	if clin and clin != "GX": 
		if "High" in clin: clin=4
		elif "1" in clin: clin=1
		elif "2" in clin: clin=2
		elif "3" in clin: clin=3
		elif "4" in clin: clin=4
		else: clin = False

		if clin:
			data[case]['histologic_grade'] = clin

	clin = load_data_from_path(dd, "tcga", "patient", 'stage_event', 'pathologic', 'text')
	
	if clin: 
		clin = clin.count('I')
		data[case]['pathologic_grade'] = clin

	
	clin = load_data_from_path(dd, "tcga", "admin:admin", 'disease_code', 'text')
	if clin:
		if clin not in disease_lookup: disease_lookup.append(clin)
		data[case]['disease'] = disease_lookup.index(clin)
	
	if not all(key in data[case] for key in keys[:-2]):
		continue

	if "mirna_expression_file" in data[case]:
		gene_data_file = data[case]['mirna_expression_file']
		df = pd.read_csv(gene_data_file, sep='\t')
		data[case]['mirna_data'] = np.array(df["reads_per_million_miRNA_mapped"])/10000.0

	if "gene_expression_file" in data[case]:
		gene_data_file = data[case]['gene_expression_file']
		df = pd.read_csv(gene_data_file, sep='\t', header=None)
		df = df[~df[0].str.contains("__")]
		data[case]['gene_data'] = np.array(df[1])



X = []
Y_recurrence = []
Y_vital_status = []
cases = []
projects = []

counts = [0] * 7
for case in data.keys():
	if not all(key in data[case] for key in keys):
		continue

	mirna_arr = data[case]['mirna_data']
	gene_arr = data[case]['gene_data']
	data_arr = np.zeros(len(gene_arr) + len(mirna_arr) + 4)
	print (len(data_arr))

	L, M = len(gene_arr), len(gene_arr) + len(mirna_arr)
	data_arr[:L] = gene_arr
	data_arr[L:M] = mirna_arr
	data_arr[M:] = np.array([data[case]['gender'], data[case]['race'], data[case]['age'], data[case]['disease']])
	project = data[case]['project']
	
	X.append(data_arr)
	Y_recurrence.append(data[case]['recurrence'])
	Y_vital_status.append(data[case]['vital_status'])

	cases.append(case)
	projects.append(project)

X = np.array(X)
Y_recurrence = np.array(Y_recurrence)
Y_vital_status = np.array(Y_vital_status)

np.savez_compressed("processed/processed_data.npz", X=X, Y_recurrence=Y_recurrence, Y_vital_status=Y_vital_status, cases=cases, projects=projects,
	race_lookup=race_lookup, drug_lookup=drug_lookup, disease_lookup=disease_lookup)


