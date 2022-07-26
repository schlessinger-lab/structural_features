STRUCTURAL FEATURES

This method generates protien structural features from genomic or proteomic data

HOW TO SET UP:

- Download the generate_structural_features.py file
- Download the database folder and unzip in the same directory as the generate_structural_features.py file
- Ensure all dependancies listed below are met

HOW TO USE:

- Make a directory in the same directory as generate_structural_features.py filled with csv files containing line separated gene or protien names, expression values
- Make an output directory in the directory as generate_structural_features.py
- In the command line run: python generate_structural_features.py input_directory_name output_directory_name
- Two additional arguments can be listed at the end of the command: True or False for the use of weights when averageing the structural features and the name of the non default background you want to use (now that you are an expert in using structural features see below for how to generate this background).

ADDITIONAL INSTRUCTIONS (now you want to get fancy):

- Creating and using your own background:

* Ensure the folder with the outputs from structural features is in the same directory as make_background.py and that the file names within have not been changed
* In the command line run: python make_background.py name_of_input_folder name_of_the_file_id_for_structural_feature_outputs name_of_output_directory
* The name_of_the_file_id_for_strucural_feature_outputs can be found by entering the directory with the output from structural features and removing the 'average_' and '.csv' from the file name
* The name_of_output_directory is the name that you will as the background_folder_name in structural features

- Updating databases with newer versions:

* If there is an update to any of the databases included in structural features do the following
* Download the database of interest
* Delete the old version of the database in the databases folder and replace it with the newly downloaded version
* Delete the files in the precounted human genome folder
* Run the update_databases.py file by typing python update_databases.py in the command line

DEPENDANCIES:
- python 3
- python packages:
* sys
* os
* pandas
* glob
* sqlite3
* math
* scipy
* statsmodels
* numpy

