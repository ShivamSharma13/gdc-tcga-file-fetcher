import argparse
import requests
import json

data_endpt = "https://api.gdc.cancer.gov/slicing/view/"

def download_file_from_gdc_api(file_uuid, token_file_path):
	file_endpt = data_endpt + file_uuid

	with open(token_file_path,"r") as token:
		token_string = str(token.read().strip())

	response = requests.post(data_endpt,
						headers = {
							"Content-Type": "application/json",
							"X-Auth-Token": token_string
							})

	file_name = file_uuid + ".bam"

	with open(file_name, "wb") as output_file:
		output_file.write(response.content)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-f", "--file-uuid", help="UUID of the file you want to download.", required=True)
	parser.add_argument("-t", "--token-file-path", help="Path to your GDC token for controlled files.", required=True)

	args = vars(parser.parse_args())

	token_file_path = args['token_file_path']
	file_uuid = args['file_uuid']

	download_file_from_gdc_api(file_uuid, token_file_path)