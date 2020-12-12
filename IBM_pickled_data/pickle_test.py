#! /usr/bin/env python

import pickle

from ClusterMapData import ClusterMapData


genera = [
    'campylobacter',
    'escherichia',
    'listeria',
    'salmonella',
    'shigella'
]

print('Pivot Domain vs. Neighbor Domain')
print('--------------------------------')
print('')

for genus in genera:
    with open(f'domain/{genus}.pickle', 'rb') as pickle_file:
        cluster_map_data = pickle.load(pickle_file)
        print(f'genus: {genus}')
        print(f'shape={cluster_map_data.data.shape}')
        print(f'row_labels={cluster_map_data.row_labels[:5]}')
        print(f'col_labels={cluster_map_data.col_labels[:5]}')
        print('')

print('Genome vs. Domain')
print('-----------------')
print('')

for genus in genera:
    with open(f'genome/{genus}.pickle', 'rb') as pickle_file:
        cluster_map_data = pickle.load(pickle_file)
        print(f'genus: {genus}')
        print(f'shape={cluster_map_data.data.shape}')
        print(f'row_labels={cluster_map_data.row_labels[:5]}')
        print(f'col_labels={cluster_map_data.col_labels[:5]}')
        print('')

