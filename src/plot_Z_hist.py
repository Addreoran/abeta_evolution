import os.path

import click
import matplotlib.pyplot as plt
import numpy as np


def generate_data(folder, files):
    position_16plus = -3
    position_total = -4
    position_to16 = -5
    result_16plus = {}
    result_total = {}
    result_to16 = {}
    possible_values = set()

    for file in files:
        result_16plus[file] = []
        result_total[file] = []
        result_to16[file] = []
        with open(os.path.join(folder + file)) as f:
            for line in f:
                if line.strip():
                    line_split = line.split("\t")
                    result_total[file].append(float(line_split[position_total]))
                    result_16plus[file].append(float(line_split[position_16plus]))
                    result_to16[file].append(float(line_split[position_to16]))
        set(result_total[file])
        possible_values = possible_values.union(set(result_total[file]))
        possible_values = possible_values.union(set(result_16plus[file]))
        possible_values = possible_values.union(set(result_to16[file]))

    possible_values = sorted(list(possible_values))
    result_total = {i: {k: (100 * list(j).count(k)) / len(j) for k in possible_values} for i, j in result_total.items()}
    result_16plus = {i: {k: (100 * list(j).count(k)) / len(j) for k in possible_values} for i, j in
                     result_16plus.items()}
    result_to16 = {i: {k: (100 * list(j).count(k)) / len(j) for k in possible_values} for i, j in result_to16.items()}
    return result_total, result_16plus, result_to16, possible_values


def plot_result(result_total, possible_values, files, save_path="1-42.svg"):
    result_total = {i: {k: (100 * list(j).count(k)) / len(j) for k in possible_values} for i, j in result_total.items()}
    barWidth = 0.12
    data = {}
    file_list = files.key()
    for name, file in files.items():
        data[name] = [result_total[file][i] for i in possible_values]
    br1 = np.arange(len(list(data.values())[0]))
    br_list = [br1]
    for f in file_list[1:]:
        br_last = br_list[-1]
        br_list.append([x + barWidth for x in br_last])
    no = 0
    colors = ['', 'r', 'g', 'b', 'c', 'y', 'orange']
    for name, data_by_group in data.items():
        plt.bar(br_list[no], data_by_group, color=colors[no], width=barWidth,
                edgecolor='grey', label=name)
        no += 1

    plt.xlabel('Animal group', fontweight='bold', fontsize=15)
    plt.ylabel('Number of Z value occur', fontweight='bold', fontsize=15)
    plt.xticks([r + barWidth for r in range(len(list(data.values())[0]))],
               possible_values)

    plt.legend()
    plt.savefig(save_path)
    plt.show()


@click.option('--final_file', default="./result/final_file_aligned.csv",
              help='Path to finale file from script "select_nonredundant_proteins.py".')
@click.option('--files', default={"Cortilago fishes": "./result/7777.csv",
                                  "Bony fishes": "./result/7898.csv",
                                  "Amphibians": "./result/8292.csv",
                                  "Birds": "./result/8782.csv",
                                  "Mammalias": "./result/40674.csv",
                                  "Reptiles": "./result/1294634.csv"
                                  },
              help='Names of files generated from count_position_z_value_table.py')
def run(final_file, files):
    result_total, result_16plus, result_to16, possible_values = generate_data(final_file, list(files.values()))

    plot_result(result_total, possible_values, files, save_path="./result/1-42.svg")
    plot_result(result_16plus, possible_values, files, save_path="./result/1-16.svg")
    plot_result(result_to16, possible_values, files, save_path="./result/16-42.svg")
