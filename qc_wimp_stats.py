import pandas as pd
import statistics
import os
import glob
import code

def read_csv(path, sample_name):
    qc_csv = path + 'basecalling_1d_' + sample_name
    wimp_csv = path + 'classification_wimp_' + sample_name
    qc_df = pd.read_csv(qc_csv)
    wimp_df = pd.read_csv(wimp_csv)
    return qc_df, wimp_df

def count_reads(wimp_df, stats_dict, species):
    readID_dict = {}
    for name in species:
        stats_dict[name]['reads'] = int(wimp_df[wimp_df['name']==name].shape[0])
        readID_dict[name] = wimp_df[wimp_df['name']==name]['readid']
    return stats_dict, readID_dict

def calculate_stats(qc_df, stats_dict, readID_dict):
    for key, value in readID_dict.items():
        count_length = []
        count_qscore = []
        #code.interact(local = dict(globals(), **locals()))
        for id in value:
            count_length.append(int(qc_df[qc_df['read_id']==id]['seqlen']))
            count_qscore.append(int(qc_df[qc_df['read_id']==id]['mean_qscore']))
        if len(count_length) == 1:
            stats_dict[key]['seqlen_mean'] = count_length[0]
            stats_dict[key]['seqlen_stdev'] = 0
        else:
            stats_dict[key]['seqlen_mean'] = statistics.mean(count_length)
            stats_dict[key]['seqlen_stdev'] = statistics.stdev(count_length)
        if len(count_length) == 1:
            stats_dict[key]['qscore_mean'] = count_qscore[0]
            stats_dict[key]['qscore_stdev'] = 0
        else:
            stats_dict[key]['qscore_mean'] = statistics.mean(count_qscore)
            stats_dict[key]['qscore_stdev'] = statistics.stdev(count_qscore)
    return stats_dict

def make_statistics_df(path, sample_name):
    qc_df, wimp_df = read_csv(path, sample_name)
    species = wimp_df['name'].unique()
    stats_dict = {name:{} for name in species}
    stats_dict, readID_dict = count_reads(wimp_df, stats_dict, species)
    stats_dict = calculate_stats(qc_df, stats_dict, readID_dict)
    stats_df = pd.DataFrame(stats_dict).transpose()
    if not os.path.exists('results'):
        os.makedirs('results')
    stats_df.to_csv('results/statistics_' + sample_name)
    #code.interact(local = dict(globals(), **locals()))

if __name__ == '__main__':
    path = 'data/WIMP/'
    files_list = glob.glob(path + 'classification_wimp_*')
    for f in files_list:
        sample_name = '_'.join(f.split('/')[-1].split('_')[2:])
        print(path, sample_name)
        make_statistics_df(path, sample_name)