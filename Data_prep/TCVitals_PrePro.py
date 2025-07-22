# Script to output each row in a tcvitals file into an individual file

# Made by Zhu (Judy) Yao on AUg 19, 2024

# Step 1: Read data from file
input_filename = "/expanse/lustre/scratch/cpruett/temp_project/Otis/TCVitals/raw_file/ep182023-tcvitals-arch.dat"

with open(input_filename, "r") as file:
    data = file.readlines()


# Step 2: Loop through each row and save to individual file
output_dir = '/expanse/lustre/scratch/cpruett/temp_project/Otis/TCVitals/'
for i, row in enumerate(data):
    row = row.strip() # Remove any leading/trailing whitespace or newline characters
    line_split = row.split() # Separate each column
    time = line_split[3]+line_split[4]
    output_filename = output_dir+time+'.Otis-tcvitals.dat'
    with open(output_filename, "w") as output_file:
        output_file.write(row)

print("Files have been created successfully!")
