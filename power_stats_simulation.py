#!/usr/bin/env python
# coding: utf-8

# In[5]:
import matplotlib

matplotlib.use('Agg')
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import os
from collections import deque
import matplotlib.pyplot as plt
import math
import seaborn as sns
import argparse
import warnings
import sys


def my_get_group(g, key):
    if key in g.groups: return True
    return False


def produce_power_samples(transcript_isoform_max, total_coverage, n_samples, output):
    for num_of_transcripts in range(1, transcript_isoform_max):
        d = deque()

        for sqrt_coverage in range(2, math.ceil(np.sqrt(total_coverage)) + 1):
            coverage = sqrt_coverage ** 2

            if ((coverage % 2) == 1):
                coverage = coverage + 1
            coverage_eR = int(coverage / 2)
            coverage_eA = coverage - coverage_eR

            delta = np.log2(coverage_eA / coverage_eR)

            # produce N samples
            for samples in range(1, n_samples):

                p_i = np.random.random(num_of_transcripts)
                p_i /= p_i.sum()

                q_i = np.random.random(num_of_transcripts)
                q_i /= q_i.sum()

                obs_prob = [p_i, q_i]
                chi2_prob, p_prob, dof_prob, expected_prob = chi2_contingency(obs_prob, correction=False)
                cohen = round(np.sqrt(chi2_prob / transcript_isoform_max), 1)

                e_Ri = np.array(np.random.multinomial(coverage_eR, p_i, 1), dtype=float)
                e_Ai = np.array(np.random.multinomial(coverage_eA, q_i, 1), dtype=float)

                for j in range(0, len(e_Ri[0])):
                    if (e_Ri[0][j] == 0 and e_Ai[0][j] == 0):
                        e_Ri[0][j] = sys.float_info.epsilon
                        e_Ai[0][j] = sys.float_info.epsilon

                obs = [e_Ri, e_Ai]

                chi2, p, dof, expected = chi2_contingency(obs)

                d.append([coverage, coverage_eR, coverage_eA, cohen, p, delta])

            input_stats = pd.DataFrame({'total_count': [item[0] for item in d], 'ref_coverage': [item[1] for item in d],
                                        'alt_coverage': [item[2] for item in d], 'cohen_w': [item[3] for item in d],
                                        'p_value': [item[4] for item in d], 'delta_log2': [item[5] for item in d]})

            ###########################################################################################################################

            subset_input_stats = input_stats.groupby(['cohen_w', 'total_count'])

            dist_values = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1]  # input_stats.cohen_w.unique()
            coverage_values = input_stats.total_count.unique()

            data = deque()

            for dist in dist_values:
                for count in coverage_values:

                    if (my_get_group(subset_input_stats, (dist, count)) == False):
                        data.append([dist, count, np.nan])
                    else:
                        if (len(subset_input_stats.get_group((dist, count))) <= 1):
                            print('warning: not enough samples for ' + ' effect size ' + str(
                                dist) + ' and total coverage ' + str(count) + ' with ' + str(
                                num_of_transcripts) + ' transcript isoforms')
                            data.append([dist, count, np.nan])
                        else:
                            rejected_null_hypothesis = subset_input_stats.get_group((dist, count))[
                                (subset_input_stats.get_group((dist, count))['p_value'] < 0.01)]
                            power = len(rejected_null_hypothesis) / len(subset_input_stats.get_group((dist, count)))
                            data.append([dist, count, power])

            power_df = pd.DataFrame(
                {'cohen_w': [item[0] for item in data], 'total_coverage_gene': [item[1] for item in data],
                 'power': [item[2] for item in data]})

            power_df.to_csv(os.path.expanduser(output + '/transcript_isoform_' + str(num_of_transcripts) + '.csv'),
                            index=False)


def power_plot(power_stats, effect_size, output):
    power_combined = pd.read_csv(power_stats)
    power_combined_subset = power_combined[(power_combined['cohen_w'] == effect_size)]

    power_combined_subset_power = power_combined_subset.groupby(['transcript_isoform', 'total_coverage_gene'])
    trans_values = power_combined_subset.transcript_isoform.unique()
    trans_values = trans_values[(trans_values > 1)]
    coverage_values = power_combined_subset.total_coverage_gene.unique()
    data = deque()

    for trans in trans_values:
        for count in coverage_values:

            if (my_get_group(power_combined_subset_power, (trans, count)) == False):
                data.append([trans, count, np.nan])
            else:
                power = power_combined_subset_power.get_group((trans, count))['power'].values[0]
                data.append([trans, count, power])

    power_df_effectSize = pd.DataFrame(
        {'transcript isoforms': [item[0] for item in data], 'total_coverage_gene': [item[1] for item in data],
         'power': [item[2] for item in data]})

    power_df_effectSize["total_coverage_gene"] = pd.Categorical(power_df_effectSize["total_coverage_gene"],
                                                                power_df_effectSize.total_coverage_gene.unique())
    data_matrix = power_df_effectSize.pivot("total_coverage_gene", "transcript isoforms", "power")

    fig = plt.figure(figsize=(8, 8))
    sns.set(font_scale=1.4)
    r = sns.heatmap(data_matrix, cmap='BuPu', cbar=False)
    r.set_xticklabels(r.get_xticklabels(), rotation=270, horizontalalignment='center')

    plt.savefig(os.path.expanduser(output + '/power_' + str(effect_size) + '.pdf'), format='pdf', dpi=3000)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    # REQUIRED
    parser.add_argument("--output_path", "--o", required=True, help="Output dir")
    parser.add_argument("--isoform_max", type=int, required=True, help="max no. of transcript isoforms")
    parser.add_argument("--effect_size", type=float, required=True, help="effect size")
    parser.add_argument("--max_coverage", type=int, required=True, choices=range(2, 1024), help="max coverage <= 1024")
    parser.add_argument("--iterations", type=int, required=False, default=1000, help="number of simulations")

    args = parser.parse_args()

    iso_counts = args.isoform_max
    output_path = args.output_path
    effect_size = args.effect_size
    max_coverage = args.max_coverage
    num_samples = args.iterations

    warnings.filterwarnings("ignore")
    script_path = sys.path[0]
    print(script_path)

    bashCommand = 'mkdir -p ' + output_path
    os.system(bashCommand)

    bashCommand = 'mkdir ' + output_path + "/temp"
    os.system(bashCommand)

    bashCommand = 'mkdir ' + output_path + "/temp/sample_files/"
    os.system(bashCommand)

    bashCommand = 'mkdir ' + output_path + "/plot/"
    os.system(bashCommand)

    print('produce samples for simuation...')
    produce_power_samples(iso_counts + 1, max_coverage + 1, num_samples, output_path + "/temp/sample_files/")

    combined_power = pd.DataFrame()

    print('combining results...')
    for count in range(1, iso_counts + 1):
        input_stats = pd.read_csv(output_path + '/temp/sample_files/transcript_isoform_' + str(count) + '.csv')
        input_stats['transcript_isoform'] = count

        combined_power = combined_power.append(input_stats)

        combined_power.to_csv(os.path.expanduser(output_path + '/power_stats_simulation.csv'), index=False)

    print('saving plot...')
    power_plot(output_path + '/power_stats_simulation.csv', effect_size, output_path + "/plot/")