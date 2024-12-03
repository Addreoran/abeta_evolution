def load_data(file):
    sequences = []
    with open(file) as f:
        for l in f:
            if l:
                # print(l.split(";")[2])
                seq = l.split(";")[1]
                sequences.append(seq)
    return sequences


def prepare_data(sequences, aa_place):
    results = {i: {j: 0 for j in
                   ["R", "H", "K", "D", "E", "S", "T", "N",
                    "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]}
               for i in aa_place}
    print(sequences)
    for seq in sequences:
        for aa in aa_place:
            print(seq)
            if seq[aa] != '-':
                results[aa][seq[aa]] += 1
    for aa in aa_place:
        results[aa] = {i: j for i, j in results[aa].items()}
    return results


import matplotlib.pyplot as plt
import numpy as np


def make_csv(data, file, sequences):
    labels = ["R", "H", "K", "D", "E", "S", "T", "N",
              "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]
    list_of_pos = data.keys()
    with open(file, "w") as f:
        f.write(f"aminoacid;{';'.join([str(i) for i in list_of_pos])}\n")
        for l in labels:
            f.write(f"{l};")
            for pos in list_of_pos:
                f.write(f"{round(data[pos][l], 2)} ({round(100 * data[pos][l] / len(sequences), 2)}%);")
            f.write("\n")


def make_plot(data, file, width=0.35, ylabel="", xlabel=""):
    labels = ["R", "H", "K", "D", "E", "S", "T", "N",
              "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]

    x = np.arange(len(labels))  # the label locations

    fig, ax = plt.subplots()
    math = ["x - 2* width / 3", "x - width / 3", "x + width / 3"]
    nr = 0
    for aa_place, values_aa in data.items():
        places = []
        for i in labels:
            places.append(values_aa[i])
        rects1 = ax.bar(eval(math[nr]), places, width, label=aa_place)
        nr += 1
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(ylabel)
    ax.set_title(xlabel)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()

    plt.show()


def run():
    aa_place = [5, 12, 13]
    data1 = load_data("../data/final_results/final_file.csv")
    data = prepare_data(data1, aa_place)
    make_csv(data, "../data/final_results/positions_mutations.csv", data1)


run()
