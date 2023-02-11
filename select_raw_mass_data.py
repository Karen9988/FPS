#!/usr/bin/python
import os
import sys
import getopt
import pandas as pd

input_file = ''
output_file = 'screened_mass_data.xlsx'
dbe = 12
mass = 400.0


def read_xlsx(file):
    raw_data_frame = pd.read_excel(file, sheet_name=0,  # assume the first sheet is raw data
                                   header=0)  # first row is header
    return raw_data_frame


def main(argv):
    """
    input file is a XLSX file that producted by Software 'Agilent MassHunter Qualitative Analysis'
    """
    global input_file
    global output_file
    global dbe
    global mass
    opts, args = getopt.getopt(argv, "f:d:m:o:", ["input_file=", "dbe=", "mass=", "output_file"])
    for opt, arg in opts:
        if opt in ('-f', '--input_file'):
            input_file = arg
        elif opt in ('-d', '--dbe'):
            dbe = int(arg)
        elif opt in ('-m', '--mass'):
            mass = int(arg)
        elif opt in ('-o', '--output_file'):
            output_file = arg

    if len(input_file) == 0:
        print('Usage: python select_raw_mass_data.py -f xxx.xlsx')
        sys.exit(-1)
    if not os.path.isfile(input_file):
        print('input file not exists, please type in correct path')
        sys.exit(-2)
    # start read file
    data_frame = read_xlsx(input_file)
    result = pd.DataFrame()
    for row in data_frame.iterrows():
        if int(row[1]['DBE']) >= dbe and float(row[1]['Mass']) >= mass:
            # print(row[1])
            result = result.append(row[1])
    result.to_excel(output_file)

if __name__ == '__main__':
    main(sys.argv[1:])
