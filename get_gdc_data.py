import requests
import json
import os

'''
This script is highly inefficient becuase GDC doesn't return full information when it comes to files data present in cases.
'''

######################################################################################
######################Instructions for running this script############################
'''
Hello.

This script gathers data from GDC (TCGA) via their API and it's companion script 
"gather_gdc_uuid_links.py" can create a UUID file per case and file that can be used to 
download raw data files.

Try to be on an educational network such as eduroam. Although I'm not certain, but I've 
seen significant delays in API response on private networks; AT&T in my case.

Step 1: Set the data_path variable to the location of your choice. It's immediately after 
	this comment section.

Step 2: Put in primary sites in primary_sites variable. Some exaples include: Kidney, 
	Ovary, Brain etc.

Step 3: In order to get complete intersection of files, put if_files as FALSE. Run this
	script once.

Step 4: Then switch back the if_files flag to TRUE and run this script once again. This
	should create some json files in your data directory. Please check.

Step 5: Run the "gather_gdc_uuid_links.py" script which will create some links file.

Note: These link files have the cases which were sequenced for both normal and tumor
	samples. 

Hope this helps.

'''

######################################################################################

data_path = 'data/'
strategies = ['RNA-Seq']
primary_sites=['Lung']
#primary_sites=['Colon']

sampe_id_to_type = {'01': 'Primary Solid Tumor', '11': 'Solid Tissue Normal', '10': 'Blood Derived Normal'}


def get_cases_from_gdc(endpt, data_size=2, json_file_name='miRNA.json', experimental_strategy_filter=None, primary_site_filter=None):
	# The 'fields' parameter is passed as a comma-separated string of single names.
	
	fields = [
		"sample_ids",
		"samples.sample_type_id",
		"files.experimental_strategy",
		"files.file_id",
		"files.data_format",
		"workflow_type"
		]

	fields = ','.join(fields)

	filters = {
		"op": "and",
		"content":[
				{
				"op": "in",
				"content":{
					"field": "cases.project.primary_site",
					"value": [primary_site_filter],
					}
				},
				{
					"op": "in",
					"content":{
						"field": "files.experimental_strategy",
						"value": [experimental_strategy_filter],
				}
			}
		]
	}

	params = {
		"filters": filters,
		"fields": fields,
		"format": "TSV",
		"size": data_size,
		}

	response = requests.post(endpt, headers = {"Content-Type": "application/json"}, json = params)

	if not response.status_code == 200:
		print("Internet Error, perhaps a query error?")
		return 0
	content = response.content.decode("utf-8")
	
	data = process_tab_seperated_data(content, json_file_name)
	
	if data is None:
		return 0

	return 1


def get_files_from_gdc(endpt, data_size=1, json_file_name='miRNA.json', experimental_strategy_filter=None, primary_site_filter=None, file_extension="BAM"):
	fields = [
		"cases.case_id",
		"cases.samples.sample_type_id",
		"data_format",
		"experimental_strategy",
		"analysis.workflow_type",
		"file_name"
		]

	fields = ','.join(fields)

	#Creating varibale to send into later function so that directory name varibales are consistent between cases and files.
	primary_site = primary_site_filter
	experimental_strategy = experimental_strategy_filter

	#Silly Solution.
	#Prostrate gland doesn't return but Prostrate returns something.
	#Splitting by " " and taking only the first word as the filter string.
	#Replacing the filter string.
	primary_site_filter = primary_site_filter.split(" ")[0]


	filters = {
		"op": "and",
		"content":[
				{
					"op": "in",
					"content":{
						"field": "cases.project.primary_site",
						"value": [primary_site],
					}
				},
				{
					"op": "in",
					"content":{
						"field": "files.experimental_strategy",
						"value": [experimental_strategy_filter],
					}
				},
				{
					"op": "in",
					"content":{
						"field": "files.data_format",
						"value": [file_extension]
					}
				}
		]
	}

	params = {
		"filters": filters,
		"fields": fields,
		"format": "TSV",
		"size": data_size,
		}

	response = requests.post(endpt, headers = {"Content-Type": "application/json"}, json = params)

	if not response.status_code == 200:
		print("Internet Error, perhaps a query error?")
		return 0
	content = response.content.decode("utf-8")

	data = add_files_data_information_to_json_files(content, primary_site, experimental_strategy)

	if data is None:
		return 0

	return 1


def organize_files_data(column_headers):
	number_of_files = 0
	for column_header in column_headers:
		#print(column_header)
		if 'file_id' in column_header:
			number_of_files+=1

	val_to_index = {
					'number_of_files': number_of_files,
					'file_ids_index': [column_headers.index('files.'+str(i)+'.file_id') for i in range(number_of_files)],
					'file_experimental_strategies_index': [column_headers.index('files.'+str(i)+'.experimental_strategy') 
													if ('files.'+str(i)+'.experimental_strategy') in column_headers 
													else None
													for i in range(number_of_files)],
					'file_data_format_index': [column_headers.index('files.'+str(i)+'.data_format') for i in range(number_of_files)],
					}
	return val_to_index


def organize_samples_data(column_headers):
	number_of_samples = 0
	for column_header in column_headers:
		if 'sample_ids' in column_header:
			number_of_samples+=1	

	val_to_index = {
					'number_of_samples': number_of_samples,
					'sample_ids_index' : [column_headers.index('sample_ids.'+str(i)) for i in range(number_of_samples)],
					'sample_type_ids_index' : [column_headers.index('samples.'+str(i)+'.sample_type_id') for i in range(number_of_samples)],
					}

	return val_to_index


def process_tab_seperated_data(content, file_name):
	'''
	content = --\t--\t--\t--\n
	The content is a tab seperated content.
	'''

	#Obtaining rows.
	content_rows = content.split('\n')
	content_rows = list(filter(None, content_rows))
	refined_content = [i.replace('\r', '').split('\t') for i in content_rows]
	
	#Getting the column header.
	column_headers = refined_content[0]
	
	#Column data
	refined_content = refined_content[1:]

	
	samples_val_to_index = organize_samples_data(column_headers)

	files_val_to_index = organize_files_data(column_headers)
	#print(val_to_index)
	#exit()

	data = []

	try:
		case_id_index = column_headers.index('id')
	except ValueError:
		return None


	for idx, row in enumerate(refined_content):
		row_in_dict = {}
		row_in_dict['case_uid'] = row[case_id_index]
		row_in_dict['sample_data'] = {'ids_present': [], 'id_types': []}
		row_in_dict['files_data'] = {'ids_present': [], 'experimental_strategies': [], 'data_format': []}

		for sample_index in samples_val_to_index['sample_ids_index']:
			if row[sample_index] is not None:
				row_in_dict['sample_data']['ids_present'].append(row[sample_index])
		
		for sample_type_index in samples_val_to_index['sample_type_ids_index']:
			if row[sample_type_index] is not None:
				row_in_dict['sample_data']['id_types'].append(row[sample_type_index])

		for file_data_format_index in files_val_to_index['file_data_format_index']:
			if row[file_data_format_index] is not None:
				row_in_dict['files_data']['data_format'].append(row[file_data_format_index])

		for file_index in files_val_to_index['file_ids_index']:
			if row[file_index] is not None:
				row_in_dict['files_data']['ids_present'].append(row[file_index])		

		for file_experimental_strategy_index in files_val_to_index['file_experimental_strategies_index']:
			if file_experimental_strategy_index is None:
				row_in_dict['files_data']['experimental_strategies'].append(None)
			else:
				row_in_dict['files_data']['experimental_strategies'].append(row[file_experimental_strategy_index])


		data.append(row_in_dict)
	

	#print(len(data))
	#print(data)
	j = json.dumps(data)

	with open(data_path + file_name, 'w') as f:
		f.write(j)

	return data


def process_json_data(txt_file_path, primary_site, strategy):
	try:
		with open(txt_file_path, 'r') as f:
			cases_data = json.load(f)
	except FileNotFoundError:
		print("No file found for {} using {}.".format(primary_site, strategy))
		return 0
	#print(len(cases_data))
	uuid_set = set()

	for case in cases_data:
		uuid_set.add(case['case_uid'])

	print('For: ' + primary_site + ' using ' + strategy)

	#Performing checks that no repeats were introduced in the JSON data files.
	if len(uuid_set) == len(cases_data):
		print(str(len(cases_data)) + " cases found.")
		print("Ok, bye.")
	else:
		print("Some repeats.")
	return 1


def add_files_data_information_to_json_files(content, primary_site, experimental_strategy):
	content_rows = content.split('\n')
	content_rows = list(filter(None, content_rows))
	refined_content = [i.replace('\r', '').split('\t') for i in content_rows]
	
	#Getting the column header.
	column_headers = refined_content[0]
	
	#Columns data
	refined_content = refined_content[1:]
	print("Number of files received: {}".format(len(refined_content)))
	for idx, column_header in enumerate(column_headers):
		if 'sample_type_id' in column_header:
			sample_type_id_index = idx
		if 'case_id' in column_header:
			case_id_index = idx
		if 'id' == column_header:
			file_id_index = idx

	json_file_name = primary_site.replace(' ', '_') + '/' + strategy + '.json'
	
	#Check if data is present for a specific experimental strategy that the loop is considering.
	try:
		with open(data_path + json_file_name, 'r') as f:
			stored_cases_data = json.load(f)
	except FileNotFoundError:
		return 0

	#Iterating over rows of files data.
	#Looks like:
	#[['BAM', '3f3128c9-1072-443c-b6ae-8b63166fa9a4', '4c8a6aee-60f7-43bc-9ff2-71f32e9836cc', '01'], 
	#['BAM', '6394e00a-1f6f-4ea2-b885-373bbb59def9', 'b663d73f-727f-487c-8bd9-fe1905ece070', '01']]
	for row in refined_content:
		#Get the case ID from the row.
		file_case_id = row[case_id_index]

		#Get the file ID from the row.
		file_id = row[file_id_index]

		#Get the sample type from the row.
		file_sample_type_id = row[sample_type_id_index]

		#Iterate over cases data to find a particular case.
		#Ids should match.

		for stored_case_idx, stored_case in enumerate(stored_cases_data):
			#Checking the IDs
			if file_case_id == stored_case['case_uid']:
				#Found the correct case.

				#Initializing an empty list to match and hold sample type data if not already present.
				try:
					new_column_for_stored_case = stored_case['files_data']['sample_type_ids']
				except KeyError:
					new_column_for_stored_case = ["" for _ in range(len(stored_case['files_data']['ids_present']))]

				for idx, stored_case_file_id in enumerate(stored_case['files_data']['ids_present']):
					if stored_case_file_id == file_id:
						#Adding the sample ID at the same index at which the file found is present in the stored data.
						new_column_for_stored_case[idx] = file_sample_type_id
				stored_cases_data[stored_case_idx]['files_data']['sample_type_ids'] = new_column_for_stored_case

	json_data = json.dumps(stored_cases_data)
	#print(json_data)
	with open(data_path + json_file_name, 'w') as f:
		f.write(json_data)

	return 1



if __name__ == '__main__':
	if_files = False
	#if_files = True

	#API endpoints.
	cases_endpt = 'https://api.gdc.cancer.gov/legacy/cases'
	files_endpt = 'https://api.gdc.cancer.gov/legacy/files'
	
	#Switch this bool if you want to replace the existing files.
	replace = True

	for primary_site in primary_sites:
		if not primary_site.replace(' ', '_') in os.listdir(data_path):
			os.mkdir(data_path+primary_site.replace(' ', '_'))
		for strategy in strategies:
			#Messed up dependency, it's path instead of file name.
			json_file_name = primary_site.replace(' ', '_') + '/' + strategy + '.json'

			if not replace:
				if strategy+'.json' in os.listdir(data_path + primary_site.replace(' ', '_')):
					print("Found " + strategy+'.json' + ' skipping...')
					continue

			if if_files:
				status = get_files_from_gdc(files_endpt, 
								data_size=10000000, 
								json_file_name=json_file_name, 
								experimental_strategy_filter=strategy, 
								primary_site_filter=primary_site,
								file_extension='BAM')
			else:						
				status = get_cases_from_gdc(cases_endpt, 
									data_size=1000000, 
									json_file_name=json_file_name, 
									experimental_strategy_filter=strategy, 
									primary_site_filter=primary_site)
			if not status:
				continue
			process_json_data(data_path + json_file_name, primary_site, strategy)

	print("All done. Bye.")













