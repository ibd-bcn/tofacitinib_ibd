import pandas as pd
import csv
import sys

# Increase the maximum field size allowed
csv.field_size_limit(sys.maxsize)


list = ["todas_w8_NR","todas_w8_R","todas_w0_R","todas_w0_NR"]

for z in range(len(list)):
    value = list[z]
    
    # Adjust the path and filename as necessary
    input_file = str('/cytosig_mine/results/'+value+'.xlsx.Zscore')  # Update this path
    output_file = str('/cytosig_mine/results/'+value+'_converted.csv')  # Desired output path

    # Use sep='\t' for tab-separated values and index_col=0 to use the first column as row names
    df = pd.read_csv(input_file, sep='\t', engine='python', index_col=0)

    # Write the data to an Excel file, retaining row names
    df.to_csv(output_file)
