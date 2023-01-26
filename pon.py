# Rough toolkit with functions for loading indexcov/mosdepth.tar.gz file data
#   and lines for prepping a panel of normals

# Libraries
import os
from glob import glob
import tarfile
import gzip
import tempfile
import argparse
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial

# TODO: TESTING VALUES =========================================================
import sys
sys.argv = ['cnv.py', '.']

# PARSE ARGUMENTS ==============================================================
parser = argparse.ArgumentParser(description = 'Call CNV in indexcov and mosdepth output data.')
parser.add_argument('data_dir',
    help = 'Input data directory.')
args = parser.parse_args()
data_set = os.path.basename(os.path.normpath(args.data_dir))
print('Processing: ' + data_set)

# ANALYSIS FUNCTIONS ===========================================================
# Load coverage data
#   Input indexcov or mosdepth. Creates temp directory to untar and gunzip
#   output files while specifically extracting the BED output file. Returns raw
#   file data.
#   - Input:
#       - args: arguments consisting of data_dir, data directory
#       - type: part string of either 'indexcov' or 'mosdepth' to choose which data to load
#       - save_name: boolean whether to insert 'dataset' column
#   - Output:
#       - cov_data: coverage data in BED format (columns chromosome, start, end, cov)
def load(data_dir, type, save_name = True):
    "Load coverage data for a single directory."
    if type in 'indexcov':
        file_type = '/indexcov.tar.gz'
    elif type in 'mosdepth':
        file_type = '/mosdepth.tar.gz'
    else:
        raise Exception("No input file type selected: mosdepth or covexcov.")
    data_set = os.path.basename(os.path.normpath(data_dir))
    with tempfile.TemporaryDirectory() as temp_dir:
        cov_tar = tarfile.open(data_dir + file_type)
        for member in cov_tar.getmembers():
            if member.name.endswith('.bed.gz'):
                cov_tar.extract(member, temp_dir)
                cov_file = member.name
                with gzip.open(temp_dir + '/' + cov_file) as cov_gz:
                    cov_data = pd.read_table(
                        cov_gz,
                        header = 0,
                        names = ['chromosome', 'start', 'end', 'depth']
                    )
                if save_name:
                    cov_data.insert(0, 'dataset', data_set)
    return(cov_data)

# Load multiple files
#   Input parent directory to load multiple files. Calls load() on list and
#   combines to one output dataframe.
#   - Inputs:
#       - parent_dir: parent directory of data files
#       - type: part string of either 'indexcov' or 'mosdepth'
#       - save_name: boolean whether to insert 'dataset' column
#   - Outputs:
#       - all_data: dataframe of all coverage data
def load_directory(parent_dir, type, save_name = True):
    "Load coverage data from all samples in a parent directory."
    folder_list = []
    [folder_list.append(folder) for folder in os.listdir(parent_dir) if 'NWD' in folder]
    with mp.Pool(mp.cpu_count()) as pool:
        run_load = partial(load, type = type, save_name = save_name)
        set_data = pool.map(run_load, folder_list)
    if save_name:
        all_data = pd.concat(set_data).sort_values(['dataset', 'chromosome', 'start']).reset_index(drop = True)
    elif not save_name:
        all_data = pd.concat(set_data).sort_values(['chromosome', 'start']).reset_index(drop = True)
    return(all_data)

# TOOLKIT FUNCTIONS ============================================================
def sample_dir(dir, n):
    "Samples directory for n NWD* folders with a fixed seed."
    folder_list = []
    [folder_list.append(folder) for folder in os.listdir(dir) if 'NWD' in folder]
    random.seed(1)
    return(random.sample(folder_list, n))

# TODO: Toolkit function to generate panel of normals
# 
def make_pon(args):
    return('PON!')

# SCRIPT =======================================================================

# Generate indexcov PON
# Load indexcov data
path_jzterra = os.environ.get('pathjzterra')
sample_dir = path_jzterra + '/outputs/pharmu/sample/'
ind_data = load_directory(sample_dir, 'indexcov', True)
#ind_data.to_csv('ind_data.csv', index = False)

# Group data into windows, flatten dataframe, save median scaled depth as PON
# TODO: Most likely does not need additional scaling since indexcov scales results
ind_win = ind_data.groupby(['chromosome', 'start', 'end'])
ind_win_stats = ind_win[['depth']].agg(['median'])
ind_win_pon = ind_win_stats.reset_index()
ind_win_pon.columns = ind_win_pon.columns.get_level_values(0)
ind_win_pon.to_csv('indexcov_pon.csv', index = False)
test = pd.read_csv('indexcov_pon.csv')

# Generate mosdepth PON
# Aggregate window data and find median to scale off of
mos_data = load_directory(sample_dir, 'mosdepth', True)
#mos_data.to_csv('mos_data.csv', index = False)

# TODO: Generate mosdepth PON
data_sets = mos_data.dataset.unique()
mos_sets = []
for dset in data_sets:
    temp_set = mos_data[mos_data.dataset == dset].reset_index(drop = True)
    temp_set['depth'] = temp_set['depth'] / temp_set['depth'].median()
    mos_sets.append(temp_set)

mos_sets = pd.concat(mos_sets).sort_values(['dataset', 'chromosome', 'start']).reset_index(drop = True)

# Group data into windows, flatten dataframe, save median scaled depth as PON
# TODO: Scale samples individually first? Then incorporate into PON? Or scale after?
mos_win = mos_sets.groupby(['chromosome', 'start', 'end'])
mos_win_stats = mos_win[['depth']].agg(['median'])
mos_win_pon = mos_win_stats.reset_index()
mos_win_pon.columns = mos_win_pon.columns.get_level_values(0)
mos_win_pon.to_csv('mosdepth_pon.csv', index = False)
test = pd.read_csv('mosdepth_pon.csv')
