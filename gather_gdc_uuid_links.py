import json
import os
#import matplotlib.pyplot as plt
import numpy as np

from get_gdc_data import strategies, primary_sites, data_path


def get_number_of_cases_by_site():
	numbers = {}
	for primary_site in primary_sites:
		numbers[primary_site] = {}
		for strategy in strategies:
			json_file_name = primary_site.replace(' ', '_') + '/' + strategy + '.json'
			if not strategy+".json" in os.listdir(data_path + primary_site.replace(' ', '_')):
				#File does not exist, skip.
				continue

			with open(data_path + json_file_name, 'r') as f:
				cases_data = json.load(f)
			numbers[primary_site][strategy] = len(cases_data)
	for k, v in numbers.items():
		print("{} {}".format(k, v))
	return numbers


def get_number_of_files(bam_filter = True):
	numbers_for_sites = {}
	for primary_site in primary_sites:
		numbers_for_strategies = {}
		for strategy in strategies:
			json_file_name = primary_site.replace(' ', '_') + '/' + strategy + '.json'
			if not strategy+".json" in os.listdir(data_path + primary_site.replace(' ', '_')):
				#File does not exist, skip.
				continue

			with open(data_path + json_file_name, 'r') as f:
				cases_data = json.load(f)	

			number_of_files_per_case = {}
			for case in cases_data:

				#Getting only the files dictionary data.
				files_data = case['files_data']
				
				#Obtaining indices from files data for a specific experimental strategy.
				exp_indices = []
				
				#print(case)
				exp_strats = files_data['experimental_strategies']
				for idx, exp_strat in enumerate(exp_strats):
					if exp_strat == strategy:
						#Only getting the .bam files.
						if not bam_filter:
							#Gathering the indices on which 
							exp_indices.append(idx)
						else:
							if files_data['data_format'][idx] == 'BAM':
								exp_indices.append(idx)
				#print(miRNA_indices)

				if len(exp_indices) not in number_of_files_per_case:
					number_of_files_per_case[len(exp_indices)] = 1
				else:
					number_of_files_per_case[len(exp_indices)] += 1

				del exp_indices
			#Making a dictionary numbers_for_stratergies which stores the number of files belonging to an 
			#experimental strategy for a specific tumor.
			numbers_for_strategies[strategy] = number_of_files_per_case
		numbers_for_sites[primary_site] = numbers_for_strategies
		del numbers_for_strategies

	#print(numbers_for_sites)
	for k, v in numbers_for_sites.items():
		print("{} {}".format(k, v))
	return numbers_for_sites


def get_sample_type_ids():
	op = {}
	for primary_site in primary_sites:
		numbers_of_samples = {}
		for strategy in strategies:
			json_file_name = primary_site.replace(' ', '_') + '/' + strategy + '.json'
			if not strategy+".json" in os.listdir(data_path + primary_site.replace(' ', '_')):
				#File does not exist, skip.
				continue
			small = {'01':0, '11':0}
			with open(data_path + json_file_name, 'r') as f:
				cases_data = json.load(f)

			for case in cases_data:
				files_data = case['files_data']
				try:
					samples_ids = files_data['sample_type_ids']
					for sample_id in samples_ids:
						try:
							small[sample_id] += 1
						except KeyError:
							pass
				except KeyError:
					pass
			numbers_of_samples[strategy] = small
		op[primary_site] = numbers_of_samples
	for k, v in op.items():
		print("{} {}".format(k, v))


def get_links_file(primary_site, strategy):
	json_data_for_output = {}
	print('{} ---- {}'.format(primary_site, strategy))

	json_file_name = primary_site.replace(' ', '_') + '/' + strategy + '.json'
	if not strategy+".json" in os.listdir(data_path + primary_site.replace(' ', '_')):
		#File does not exist, skip.
		return 0
	with open(data_path + json_file_name, 'r') as f:
		cases_data = json.load(f)

	#Looping over the cases in the json file.
	for case in cases_data:
		#Getting the case ID.
		case_id = case['case_uid']

		json_data_for_output[case_id] = {}

		#Getting the dictionary files_data.
		files_data = case['files_data']

		#Filtering based on experimental_strategey and data_format of all files stored for a specific case data.
		for idx, experimental_strategey in enumerate(files_data['experimental_strategies']):
			#Check if experimental strategy matches the strategy in focus.
			if experimental_strategey == strategy:
				#Check if the file is a BAM file.
				if files_data['data_format'][idx] == 'BAM':
					#Storing the file ID and tumor ample type of the file
					try:
						#In case the list for tumor sample type ID doesn't exist.
						json_data_for_output[case_id][files_data['ids_present'][idx]] = files_data['sample_type_ids'][idx]
					except KeyError:
						continue

	return json_data_for_output


def flush_loner_strings(file_string):
	'''
	Function made to delete the loner case IDs which made to the final file.
	We're only considering 01s and 11s in our probelem.
	'''

	rows = file_string.split('\n')
	data = [row.split('\t') for row in rows]
	data = list(filter(None, data))
	checker = {}
	for case in data:
		case_id = case[0]
		sample_type = case[-1]
		if case_id is '' or case is '' or sample_type is '':
			#Don't consider empty/trash case ids.
			continue
		if case_id not in checker:
			checker[case_id] = {'01': 0, '11': 0}
		checker[case_id][sample_type] += 1

	bad_cases = [case_id for case_id, sample_type_numbers in checker.items() if 0 in sample_type_numbers.values()]
	
	#Let's remove the bad cases.
	remove_idxs = []
	for idx, case in enumerate(data):
		case_id = case[0]
		for bad_case in bad_cases:
			if case_id == bad_case:
				remove_idxs.append(idx)

	output = ['\t'.join(row) for idx, row in enumerate(data) if idx not in remove_idxs]
	
	return '\n'.join(output)

def	filter_and_store(primary_site, strategy, json_data):
	'''
	Looks like:
	{'00e41e6a-9fe7-44f9-978b-7b05b179506a': {'62389140-264f-4011-9bad-ce812d871e19': '01', 
											'eec525cf-3c53-4767-87d6-79516a96413d': '11'},
											}
	'''
	file_string = ''
	for case_id, files_data in json_data.items():
		if len(files_data) == 1:
			pass
		elif len(files_data) > 1:
			for file_id, sample_type in files_data.items():
				if sample_type == '01' or sample_type == '11':
					file_string += '{}\t{}\t{}\n'.format(case_id, file_id, sample_type)
				else:
					break

	links_file_name = data_path + primary_site.replace(' ', '_') + '/' + strategy + '_links' + '.txt'
	
	file_string = flush_loner_strings(file_string)

	with open(links_file_name, 'w') as f:
		f.write(file_string)
	
	print("Done: {}".format(links_file_name))
	return 1


def backbone():
	'''
	Getting number of cases with experimental strategies.
	'''
	#case_numbers_by_sites = get_number_of_cases_by_site()
	'''
	for site, numbers in case_numbers_by_sites.items():
		
		create_pie_charts([val for key,val in numbers.items()],
						[key+" (Number of cases: "+str(val)+")" for key,val in numbers.items()],
						site
						)
	'''
	
	#file_number_distribution = get_number_of_files()
	
	#get_sample_type_ids()
	for primary_site in primary_sites:
		for strategy in strategies:
			json_data = get_links_file(primary_site, strategy)
			if not json_data:
				continue
			filter_and_store(primary_site, strategy, json_data)



if __name__ == "__main__":
	backbone()	

