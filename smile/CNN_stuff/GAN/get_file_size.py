import os
from gan_secrets import csv_file

def get_file_size(file_path):
    """
    Get the size of a file in bytes and convert to human-readable format
    """
    if not os.path.isfile(file_path):
        return f"File not found: {file_path}"
    
    # Get file size in bytes
    size_bytes = os.path.getsize(file_path)
    
    # Convert to human-readable format
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0 or unit == 'GB':
            break
        size_bytes /= 1024.0
    
    return f"File: {file_path}\nSize: {size_bytes:.2f} {unit}"

def get_dataframe_info(file_path):
    """
    Get information about the CSV file as a pandas DataFrame
    """
    import pandas as pd
    
    try:
        df = pd.read_csv(file_path)
        rows, cols = df.shape
        return f"DataFrame dimensions: {rows} rows × {cols} columns"
    except Exception as e:
        return f"Error reading CSV: {str(e)}"

if __name__ == "__main__":
    print(get_file_size(csv_file))
    print(get_dataframe_info(csv_file))