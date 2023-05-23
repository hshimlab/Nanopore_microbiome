'''
2023.08.02
Yunseol Park
'''

import pandas as pd
import statistics
import os
import glob

def read_csv(path, sample_name):
    '''
    Function that reads in the two files (basecalling and WIMP).

    Args:
        path (str): directory name of where the WIMP files are stored.
        sample_name (str): name of the sample
    Returns:
        qc_df (DataFrame): data frame of the basecalling file
        wimp_df (DataFrame): data frame of the WIMP file
    '''
    qc_csv = path + 'basecalling_1d_' + sample_name
    wimp_csv = path + 'classification_wimp_' + sample_name
    qc_df = pd.read_csv(qc_csv)
    wimp_df = pd.read_csv(wimp_csv)
    return qc_df, wimp_df

def count_reads(wimp_df, stats_dict, species):
    '''
    Function that counts the number of reads.

    Args:
        wimp_df (DataFrame): data frame of the WIMP file
        stats_dict (dict): nested dictionary that has species name as key and statistics as value
        species (list): list of species names
    Returns:
        stats_dict (dict): nested dictionary that has species name as key and statistics as value
        readID_dict (dict): dictionary containing the species name as key and the read ID as value
    '''
    readID_dict = {}
    for name in species:
        # Get number of reads
        stats_dict[name]['reads'] = int(wimp_df[wimp_df['name']==name].shape[0])
        # Get read ID
        readID_dict[name] = wimp_df[wimp_df['name']==name]['readid']
    return stats_dict, readID_dict

def calculate_stats(qc_df, stats_dict, readID_dict):
    '''
    Function that calculates the statistics.

    Args:
        qc_df (DataFrame): data frame of the basecalling file
        stats_dict (dict): nested dictionary that has species name as key and statistics as value
        readID_dict (dict): dictionary containing the species name as key and the read ID as value
    Returns:
        stats_dict (dict): nested dictionary that has species name as key and statistics as value
    '''
    for key, value in readID_dict.items():
        # Get the read lengths and the quality score
        count_length = []   # Length of sequence (read)
        count_qscore = []   # Quality score
        for id in value:
            count_length.append(int(qc_df[qc_df['read_id']==id]['seqlen']))
            count_qscore.append(int(qc_df[qc_df['read_id']==id]['mean_qscore']))
        # Calculate the mean of read length
        if len(count_length) == 1:
            stats_dict[key]['seqlen_mean'] = count_length[0]
            stats_dict[key]['seqlen_stdev'] = 0
        else:
            stats_dict[key]['seqlen_mean'] = statistics.mean(count_length)
            stats_dict[key]['seqlen_stdev'] = statistics.stdev(count_length)
        # Caluclate the Q score mean
        if len(count_length) == 1:
            stats_dict[key]['qscore_mean'] = count_qscore[0]
            stats_dict[key]['qscore_stdev'] = 0
        else:
            stats_dict[key]['qscore_mean'] = statistics.mean(count_qscore)
            stats_dict[key]['qscore_stdev'] = statistics.stdev(count_qscore)
    return stats_dict

def make_statistics_df(path, sample_name):
    '''
    Function that calculates statistics of WIMP data and writes to CSV.

    Args:
        path (str): directory name of where the WIMP files are stored.
        sample_name (str): name of the sample
    Returns:
        None
    '''
    # Read in the data
    qc_df, wimp_df = read_csv(path, sample_name)
    # Get all species names
    species = wimp_df['name'].unique()
    # Dictionary to save the statistics
    stats_dict = {name:{} for name in species}
    # Get number of reads and read ID
    stats_dict, readID_dict = count_reads(wimp_df, stats_dict, species)
    # Calculate statistics and convert into dataframe
    stats_dict = calculate_stats(qc_df, stats_dict, readID_dict)
    stats_df = pd.DataFrame(stats_dict).transpose()
    # Save to CSV
    if not os.path.exists('results'):
        os.makedirs('results')
    stats_df.to_csv('results/statistics_' + sample_name)

if __name__ == '__main__':
    path = 'data/WIMP/'
    files_list = glob.glob(path + 'classification_wimp_*')
    for f in files_list:
        sample_name = '_'.join(f.split('/')[-1].split('_')[2:])
        print(path, sample_name)
        make_statistics_df(path, sample_name)