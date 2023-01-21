# Author: Rediet
# Step 6 of Smile. This file calculates the minium spectral angle. 

def determine_min_sa(sa_deg):
  """Calculates the minimum spectral angle
    Variables used: 
    min_each_row[g_data_dim[2]] holds the minimum value of each row in sa_deg (comparing between shifts for each column of the data)
    min_col_num[g_data_dim[2]] holds the colmumn number of the minimum values of each row in sa_deg

  Args: 
    sa_deg: A 2D matrix (g_data_dim[2], g_total_shifts)containing SA in degrees 

  Returns: 
    min_col_num: A 1D matrix of size g_data_dim[2] containing the shifts of the best matched spectrum in sa_deg
  """
  #Create a 1d array of size g_data_dim[2] containing the smallest value in each row of sa_deg (the spectral angle in degrees)
  #Create another 1D array of size g_data_dim[2] containing the column number of the smallest value in each row
  min_each_row = np.zeros(g_data_dim[2])
  min_col_num = np.zeros(g_data_dim[2])

  for i in range(0, g_data_dim[2]):
    col = sa_deg[i]
    #seaching each row for the minimum 
    min_each_row[i] = min(col)
    #finding the first occurance of the index of min_each_row[i]
    min_col_num[i] = np.where(col == min_each_row[i])[0][0]  
  
  return min_col_num.astype(int)
