'''
Read the input file.

The input file is a tab delimited file with the first column indicating the
sample ID, the second column indicates the output that one wants to predict,
the third column indicates whether the data should be treated as training, or
testing set, and the next two columns must have the path to the read files.

Assumes PE data, but allow for SE.
'''

import csv
import os
import sys

def open_file(filename):
    '''
    Helper function to open a file. Throw an error if not able to read the file
    '''
    try:
        filename = os.path.abspath(filename)
        fin = open(filename, 'r')
        return(fin)
    except IOError:
        print("Could not open file {}".format(filename))

def read(filename):
    '''
    Helper function to read an input file.

    Tries to guess the delimiter.
    '''
    fin = open_file(filename)
    # try to guess the field delimiter
    dialect = csv.Sniffer().sniff(fin.read(2048))
    fin.seek(0)
    reader = csv.reader(fin, delimiter = dialect.delimiter, skipinitialspace = True)
    tab = [row for row in reader]
    fin.close()
    return tab

def check_paths(tab):
    pass
