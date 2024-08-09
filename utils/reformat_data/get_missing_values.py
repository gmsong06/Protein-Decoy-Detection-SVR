import csv

def save_pdb_files_to_txt(csv_filename, txt_filename):
    # Open the CSV file
    with open(csv_filename, mode='r', newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        
        # Open the text file in write mode
        with open(txt_filename, mode='w') as txt_file:
            # Loop through each row in the CSV file
            for row in csv_reader:
                # Get the value from the 'pdb_file' column
                pdb_file_value = row['pdb_file']
                # Write the value to the text file
                txt_file.write(pdb_file_value + '\n')

# Example usage
csv_filename = '4u1c/4u1c_all_contacts.csv'
txt_filename = 'missing_pdb_files.txt'
save_pdb_files_to_txt(csv_filename, txt_filename)
