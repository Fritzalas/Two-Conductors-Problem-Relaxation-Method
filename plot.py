import os
import pandas as pd

# Specify the directory where your DAT files are located
directory = "/home/kali/Ypologistiki_Fysiki_I/ge18015-FRITZALAS-KOSMAS"

# Initialize an empty DataFrame to store the data
combined_data = pd.DataFrame()

# Loop through all DAT files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".dat"):
        file_path = os.path.join(directory, filename)
        
        # Read the DAT file into a DataFrame without changing columns
        data = pd.read_csv(file_path, sep='\t', header=None)
        
        # Append the data to the combined_data DataFrame
        combined_data = pd.concat([combined_data, data])

# Reset the index of the combined data
combined_data.reset_index(drop=True, inplace=True)

# Specify the path for the new DAT file
output_dat_file = "combined_data.dat"

# Save the combined data to a new DAT file
combined_data.to_csv(output_dat_file, sep='\t', header=False, index=False)

print(f"Combined data saved to {output_dat_file}")


