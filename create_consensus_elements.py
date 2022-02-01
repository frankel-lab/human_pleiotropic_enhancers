#!/usr/bin/env python
import pandas as pd
import sys

file_name = sys.argv[1]

summits = pd.read_csv(file_name, header=None, sep='\t')

def create_consensus_elements(ranges):
    '''Create consensus elements from summit cluster file'''
    ranges = ranges.split(',')
    x = [(j.split('-')[1:]) for j in ranges]
    chromosome = [(j.split('-')[0]) for j in ranges]

    flat_list = [int(item) for sublist in x for item in sublist]
    return f'{chromosome[0]}\t{min(flat_list)}\t{max(flat_list)}\n'

x = summits[3].apply(lambda x: create_consensus_elements(x))

with open('consensus_elements.bed', 'wt') as f:
    for j in x:
        f.write(j)
