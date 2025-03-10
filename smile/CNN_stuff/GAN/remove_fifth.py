# Very simple code to get rid of the 5th column of the dataset. It says training and testing, but to allow for better training, we're removing it 
import pandas as pd
from gan_secrets import csv_file

# Read the CSV file
data = pd.read_csv(csv_file)

# Display information about the original data
print(f"Original data shape: {data.shape}")
print(f"Original columns: {data.columns.tolist()}")

# Using column index
fifth_column_name = data.columns[4]
data_modified = data.drop(columns=[fifth_column_name])

# Display information about the modified data
print(f"Modified data shape: {data_modified.shape}")
print(f"Modified columns: {data_modified.columns.tolist()}")
print(f"Removed column: {fifth_column_name}")

# Replace original data with modified data
output_file = csv_file

# Setting index to false makes it so the index is not included within the data frame (doesn't really matter if its included or not, but you have to choose one of them)
'''
index = True
, value
0, 123
1, 213
3, 543

index = False
value 
123
213
543
'''
data_modified.to_csv(output_file, index=False)
print(f"Modified data saved to {output_file}")