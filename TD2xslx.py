#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 14:33:05 2023

@author: msharp
"""

import pandas as pd

#excel_file_path = "/Users/msharp/Desktop/Apr_May_TD_statement.xlsx"
#output_path = "/Users/msharp/Desktop/Apr_May_TD_statement_split.xlsx"

excel_file_path = "/Users/msharp/Desktop/Mar_Apr_TD_statement.xlsx"
output_path = "/Users/msharp/Desktop/Mar_Apr_TD_statement_split.xlsx"

df = pd.read_excel(excel_file_path)


# Replace the 2nd space with a comma
df['PDF_Info'] = df['PDF_Info'].replace('^([^ ]* [^ ]*) ', r'\1;', regex=True)

# Replace the 3rd (formerly 4th) space with a comma
df['PDF_Info'] = df['PDF_Info'].replace('^([^ ]* [^ ]* [^ ]*) ', r'\1;', regex=True)

# Convert " $" to ",$"
df['PDF_Info'] = df['PDF_Info'].str.replace(' $', ';$', regex=False)

# Convert "-$" to ",,$"
df['PDF_Info'] = df['PDF_Info'].str.replace('-$', ';;$', regex=False)

# Split columns using ","
df_split = df['PDF_Info'].str.split(';', expand=True)

# Get rid of rows that don't start with APR or MAY
df_split = df_split[df_split[0].str.startswith(('MAR', 'APR', 'MAY'))]

# Save the updated DataFrame back to Excel if needed
df_split.to_excel(output_path, index=False)

