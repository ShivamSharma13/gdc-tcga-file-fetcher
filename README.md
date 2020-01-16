<h1>File UUID fetcher for GDC/TCGA</h1>

<h4>These two scripts together will fetch all the cases in the GDC which have Primary Solid Tumor and Solid Tissue Normal sample sequenced.</h4>
<h4>You can put in filters for type of sequencing and type of tumor.</h4>

<h4>Detailed instructions of how to use these scripts are present at the start of "get_gdc_data.py" script.</h4>


<h4>Try to be on an educational network such as eduroam. Although I'm not certain, but I've 
seen significant delays in API response on private networks; AT&T in my case.

Step 1: Set the data_path variable to the location of your choice. It's immediately after 
	this comment section.

Step 2: See what all sequencing strategies you want, and set the strategies variable 
	accordingly.

Step 3: Put in primary sites in primary_sites variable. Some exaples include: Kidney, 
	Ovary, Brain etc.

Step 4: In order to get complete intersection of files, put if_files as FALSE. Run this
	script once.

Step 5: Then switch back the if_files flag to TRUE and run this script once again. This
	should create some json files in your data directory. Please check.

Step 6: Run the "gather_gdc_uuid_links.py" script which will create some links file.

Note: These link files have the cases which were sequenced for both normal and tumor
	samples. 

Hope this helps! Thanks.</h4>
