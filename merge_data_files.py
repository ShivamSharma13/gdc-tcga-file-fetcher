import argparse
#This script will merge different sequencing technique datasets: RNA, miRNA & WXS.


def read_and_remove_repeats(data_files):
	output = {'rna': None, 'mirna': None, 'wxs': None}
	for sequencing_type, file in data_files.items():
		with open(file) as f:
			raw = f.read()
		rows = [i.split("\t") for i in raw.split("\n") if i != ""][1:]
		data = {i[0]: {'01': [], '11':[]} for i in rows}
		print("Original Data File {} Length: {}".format(sequencing_type, len(rows)))

		for entry in rows:
			case_id, barcode, sample_type, size, file_id = entry[0], entry[1], entry[2], entry[3], entry[4], 
			#Expected data struct: {case_id: {'01': [barcode, file_id, size], '11':[barcode, file_id, size]}}		
			if len(data[case_id][sample_type]) == 0:
				data[case_id][sample_type] = [barcode, file_id, size]
				continue
			else:
				if size > data[case_id][sample_type][2]:
					data[case_id][sample_type] = [barcode, file_id, size]
				else:
					#We want the largest sample data.
					continue
		output[sequencing_type] = data
		del data
		del rows

	return output


def merge_data(data, rna_file, mirna_file, wxs_file):
	rna_cases = set(data['rna'].keys())	
	mirna_cases = set(data['mirna'].keys())	
	wxs_cases = set(data['wxs'].keys())	

	#Intersection of three lists.
	common_cases = list(rna_cases & mirna_cases & wxs_cases)
	rna_intersection_entries = []
	mirna_intersection_entries = []
	wxs_intersection_entries = []
	for common_case in common_cases:
		#Expected Data Struct intersection_entries = [case_id, file_id, sample_type, size, barcode]
		##Get RNA files.
		if len(data['rna'][common_case]['01']) != 0:
			barcode, file_id, size = data['rna'][common_case]['01']
			rna_intersection_entries.append([common_case, file_id, 'Tumor', size, barcode])
		if len(data['rna'][common_case]['11']) != 0:
			barcode, file_id, size = data['rna'][common_case]['11']
			rna_intersection_entries.append([common_case, file_id, 'Normal', size, barcode])
		
		##Get miRNA files.
		if len(data['mirna'][common_case]['01']) != 0:
			barcode, file_id, size = data['mirna'][common_case]['01']
			mirna_intersection_entries.append([common_case, file_id, 'Tumor', size, barcode])
		if len(data['mirna'][common_case]['11']) != 0:
			barcode, file_id, size = data['mirna'][common_case]['11']
			mirna_intersection_entries.append([common_case, file_id, 'Normal', size, barcode])
		
		##Get WXS files.
		if len(data['wxs'][common_case]['01']) != 0:
			barcode, file_id, size = data['wxs'][common_case]['01']
			wxs_intersection_entries.append([common_case, file_id, 'Tumor', size, barcode])
		if len(data['rna'][common_case]['11']) != 0:
			barcode, file_id, size = data['wxs'][common_case]['11']
			wxs_intersection_entries.append([common_case, file_id, 'Normal', size, barcode])

	#Sort the lists by case uuids.
	rna_intersection_entries.sort(key=lambda x: x[0])
	mirna_intersection_entries.sort(key=lambda x: x[0])
	wxs_intersection_entries.sort(key=lambda x: x[0])

	#Write the output files.
	output_rna_file = rna_file.replace(".txt", ".filtered.txt")
	output_mirna_file = mirna_file.replace(".txt", ".filtered.txt")
	output_wxs_file = wxs_file.replace(".txt", ".filtered.txt")

	output_sleuth_meta_rna_file = rna_file.replace(".txt", ".sleuth_meta.txt")
	output_sleuth_meta_mirna_file = mirna_file.replace(".txt", ".sleuth_meta.txt")
	output_sleuth_meta_wxs_file = wxs_file.replace(".txt", ".sleuth_meta.txt")

	with open(output_rna_file, "w") as f:
		f.write("\n".join(["\t".join(i) for i in rna_intersection_entries]))

	with open(output_mirna_file, "w") as f:
		f.write("\n".join(["\t".join(i) for i in mirna_intersection_entries]))

	with open(output_wxs_file, "w") as f:
		f.write("\n".join(["\t".join(i) for i in wxs_intersection_entries]))

	#Sort the lists by file uuids.
	rna_intersection_entries.sort(key=lambda x: x[1])
	mirna_intersection_entries.sort(key=lambda x: x[1])
	wxs_intersection_entries.sort(key=lambda x: x[1])

	with open(output_sleuth_meta_rna_file, "w") as f2:
		f2.write("sample\tcondition\n")
		f2.write("\n".join(["\t".join([i[1], i[2]]) for i in rna_intersection_entries]))

	with open(output_sleuth_meta_mirna_file, "w") as f2:
		f2.write("sample\tcondition\n")
		f2.write("\n".join(["\t".join([i[1], i[2]]) for i in mirna_intersection_entries]))
	
	with open(output_sleuth_meta_wxs_file, "w") as f2:
		f2.write("sample\tcondition\n")
		f2.write("\n".join(["\t".join([i[1], i[2]]) for i in wxs_intersection_entries]))

	return rna_intersection_entries, mirna_intersection_entries, wxs_intersection_entries


def main():
	parser = argparse.ArgumentParser()

	#Arguments.
	parser.add_argument("-r", "--rna-file", help="Path to RNA-Seq file", required=True)
	parser.add_argument("-m", "--mirna-file", help="Path to miRNA-Seq file", required=True)
	parser.add_argument("-x", "--wxs-file", help="Path to WXS file", required=True)


	#Parse and gather whatever the user sent.
	args = vars(parser.parse_args())
	rna_file = args['rna_file']
	mirna_file = args['mirna_file']
	wxs_file = args['wxs_file']


	data = read_and_remove_repeats({'rna': rna_file, 'mirna': mirna_file, 'wxs': wxs_file})
	merge_data(data, rna_file, mirna_file, wxs_file)


if __name__ == "__main__":
	main()