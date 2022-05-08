#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import os
from collections import deque
from os import listdir
from os.path import isfile, join
import argparse
import warnings
import sys


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def f(df):
    return df[df['total_coverage'] == np.max(df.total_coverage)]


###############################################################################################################################

def apply_power_real_data(threshold, power_simulation_file, data_dir, output_dir):
    power_combined = pd.read_csv(power_simulation_file)

    mypath = (data_dir)
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    for file in files:
        if (file.startswith('._')):
            continue
        print(file)
        hap_counts_all = pd.read_csv(mypath + file, sep='\t')  # files should be tab-delimited

        hap_counts = hap_counts_all[(hap_counts_all['Alt_count'] >= 10) | (hap_counts_all['Ref_count'] >= 10)]

        hap_counts['combine'] = hap_counts['Gene'] + "_" + hap_counts['variantID']

        gene_var_ids = hap_counts['combine'].unique()

        power_combined_subset = power_combined[(power_combined['cohen_w'] == threshold)]
        coverage_values = power_combined_subset.total_coverage_gene.unique()

        data = deque()
        for gene_var in gene_var_ids:

            hap_count_subset = hap_counts[(hap_counts['combine'] == gene_var)]
            number_of_isoforms = len(hap_count_subset)
            gene_id = gene_var.split('_')[0]
            var_id = gene_var.split('_', 1)[1]

            if (number_of_isoforms == 1):
                print("gene_dropped", gene_var)
                continue

            if ((hap_count_subset['Alt_count'] == 0).all() or (hap_count_subset['Ref_count'] == 0).all()):
                print("gene_dropped", gene_var)
                continue

            total_coverage = hap_count_subset['Ref_count'].sum() + hap_count_subset['Alt_count'].sum()

            closest_coverage = closest(coverage_values, total_coverage)

            power_row = power_combined_subset.loc[
                (power_combined_subset['transcript_isoform'] == number_of_isoforms) & (
                        power_combined_subset['total_coverage_gene'] == closest_coverage), 'power']

            if (len(power_row) == 0):
                power = np.nan
            else:
                power = power_row.values[0]

            real_counts = [[hap_count_subset['Ref_count']], [hap_count_subset['Alt_count']]]

            chi2, p_value_real, dof, expected = chi2_contingency(real_counts)

            data.append([gene_id, var_id, total_coverage, power, p_value_real])

            output_power = pd.DataFrame(
                {'gene_id': [item[0] for item in data], 'variant_id': [item[1] for item in data],
                 'total_coverage': [item[2] for item in data], 'power': [item[3] for item in data],
                 'p_value': [item[4] for item in data]})

            # select the SNPwith the highest sum count in case of multiple SNP for a gene
            output_power_filtered = (output_power.groupby("gene_id").apply(f))
            output_power_filtered_indexed = output_power_filtered.set_index('gene_id')

            output_power_filtered_indexed.to_csv(output_dir + file, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # REQUIRED
    parser.add_argument("--output_dir", "--o", required=True, help="Output dir")
    parser.add_argument("--threshold", type=float, required=True, help="The effect size threshold")
    parser.add_argument("--data_dir", required=True, help="input dir")
    parser.add_argument("--power_simulation_file", required=True,
                        help="power simulation file which is the output of running power_stats_simulation")

    args = parser.parse_args()

    threshold = args.threshold
    output_dir = args.output_dir
    data_dir = args.data_dir
    power_simulation_file = args.power_simulation_file

    warnings.filterwarnings("ignore")
    script_path = sys.path[0]
    print(script_path)

    bashCommand = 'mkdir -p ' + output_dir
    os.system(bashCommand)

    print('...')
    apply_power_real_data(threshold, power_simulation_file, data_dir, output_dir)
