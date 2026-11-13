#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This convert ECMWF csv file to get aka bika coefficient to interpolate with CDO

Authors
Paolo Davini (CNR-ISAC, Nov 2023)
"""

import pandas as pd

verticals = ['L19', 'L31', 'L62', 'L91']

for vertical in verticals:
    # Load the CSV file into a Pandas DataFrame
    df = pd.read_csv(f'grids/{vertical}.csv')

    # Extract the 'a' and 'b' coefficients from the first row of the DataFrame
    a_coefficient = df.at[0, 'a [Pa]']
    b_coefficient = df.at[0, 'b']

    # Define the output file name
    output_file = f'grids/{vertical}.txt'

    # Create the formatted content
    formatted_content = f"0\t{a_coefficient:.17f}\t{b_coefficient:.17f}\n"

    for i in range(1, len(df)):
        a_value = df.at[i, 'a [Pa]']
        b_value = df.at[i, 'b']
        formatted_content += f"{i}\t{a_value:.17f}\t{b_value:.17f}\n"

    # Write the formatted content to the output file
    with open(output_file, 'w', encoding='utf8') as file:
        file.write(formatted_content)

    print(f'Output file "{output_file}" has been created.')