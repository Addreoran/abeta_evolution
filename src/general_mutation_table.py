import click

AMINOACIDS = ["R", "H", "K", "D", "E", "S", "T", "N",
              "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W", "X"]


def load_data(file):
    sequences = []
    with open(file) as f:
        for l in f:
            if l:
                seq = l.split(";")[1]
                sequences.append(seq)
    return sequences


def prepare_data(sequences, aa_place):
    results = {i: {j: 0 for j in AMINOACIDS}
               for i in aa_place}
    for seq in sequences:
        for pos in aa_place:
            if seq[pos] != '-':
                results[pos][seq[pos]] += 1
    return results


def make_csv(data, file, sequences):
    labels = AMINOACIDS
    list_of_pos = data.keys()
    with open(file, "w") as f:
        f.write(f"aminoacid;{';'.join([str(i) for i in list_of_pos])}\n")
        print(
            """\multicolumn{1}{c}{Amino acid} & \multicolumn{2}{c}{Amino acid position}\\
            \cline{2-4}
            \multicolumn{1}{c}{type} & 5 & 12 & 13\\
            \midrule
            """)
        for l in labels:
            f.write(f"{l};")
            print(f"{l}", end=" &")
            tolist = ""
            for pos in list_of_pos:
                f.write(
                    f"{round(data[pos][l], 2)};".replace("%", "\\%"))
                tolist += f"{round(data[pos][l], 2)} ({round(100 * data[pos][l] / len(sequences), 2)}%) ".replace("%",
                                                                                                                  "\\%") + "|"
            tolist = tolist.replace("|", "&", 2)
            tolist = tolist.replace("|", "")
            print(tolist + "\\\\")
            # print(''' \\\\''')
            f.write("\n")
        print(
            """
\\\\botrule
\end{tabular}}{}
\end{minipage}
\end{center}
\end{table}
            """)


@click.option('--final_file', default="./result/final_file_aligned.csv",
              help='Path to file with all TaxIds and peptides. File generated with script select_non_redundant_proteins.py')
@click.option('--aminoacids_positions', default=[5, 12, 13],
              help='Positions of aminoacids, that should be in table.')
@click.option('--result_file', default="./result/positions_mutations.csv",
              help='Path to file with result data.')
def run(aminoacids_positions, final_file, result_file):
    data1 = load_data(final_file)
    data = prepare_data(data1, aminoacids_positions)
    make_csv(data, result_file, data1)


if __name__ == "__main__":
    run()
