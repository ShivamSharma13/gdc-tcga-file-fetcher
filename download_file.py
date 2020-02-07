import argparse
import requests
import json
import os
import subprocess

data_endpt = "https://api.gdc.cancer.gov/data/"


def process_links_file(links_file):
	#Check if the file exists.
	if not os.path.exists(links_file):
		print("Links file does not exist")
		return False

	#Dict looks like: {'case_id' : {'file_id': '01', 'file_id': '11'}}
	case_to_file = {}
	
	#Read the links file.
	with open(links_file) as f:
		raw = f.read()

	rows = raw.split('\n')
	#Case ID in the first column, file UUIDs are stored in the second column, type in 3rd.

	for row in rows:
		cols = row.split('\t')
		if len(cols) < 3:
			continue
		case_id, file_id, file_type = cols[0], cols[1], cols[2]
		if case_id in case_to_file:
			case_to_file[case_id][file_id] = file_type
		else:
			case_to_file[case_id] = {file_id: file_type}

	return case_to_file


def download_file_from_gdc_api(file_uuid, token_file_path):
	file_endpt = data_endpt + file_uuid

	with open(token_file_path,"r") as token:
		token_string = str(token.read().strip())

	#print(file_endpt)

	response = requests.get(file_endpt, headers = {"Content-Type": "application/json", "X-Auth-Token": token_string})

	file_name = file_uuid + ".bam"

	with open(file_name, "wb") as output_file:
		output_file.write(response.content)


def download_file_from_gdc_client(file_id, token_file_path, output_sub_directory):
	try:
		gdc_output = subprocess.check_output(["gdc-client", "download", file_id, "-d", output_sub_directory, "-t", token_file_path])
	except subprocess.CalledProcessError:
		print("GDC client failed to run job for file UUID: {}".format(file_id))

	return

def download_manager(case_to_file, output_directory, token_file_path):
	#Make sure an output directory exists, if not the create one.
	if not os.path.exists(output_directory):
		os.mkdir(output_directory)

	for case, files in case_to_file.items():
		#Check or create a subdirectory.
		case_subdir_path = output_directory.rstrip('/') + '/' + case
		if not os.path.exists(case_subdir_path):
			os.mkdir(case_subdir_path)

		for file_id, file_type in files.items():
			#Create a subdir if required for file type.
			if not os.path.exists(case_subdir_path + '/' + file_type):
				os.mkdir(case_subdir_path + '/' + file_type)

			download_file_from_gdc_client(file_id, token_file_path, case_subdir_path + '/' + file_type)

	return True


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-f", "--links-file", help="Path to the links file you want to download.", required=True)
	parser.add_argument("-t", "--token-file-path", help="Path to your GDC token for controlled files.", required=False)
	parser.add_argument("-o", "--output-directory", help="Path to your output directory.", required=True)

	parser.add_argument("-g", "--use-gdc-client", help="Use gdc-client tool to download files.", action="store_true")

	args = vars(parser.parse_args())

	token_file_path = args['token_file_path']
	links_file = args['links_file']
	output_directory = args['output_directory']
	if_gdc = args['use_gdc_client']

	case_to_file = process_links_file(links_file)
	
	if case_to_file is False:
		exit()

	download_manager(case_to_file, output_directory, token_file_path)

