import os.path

from ete3 import NCBITaxa
import click


def read_file(path):
    result = {}
    with open(path) as f:
        for l in f:
            if l.strip():
                line = l.strip().split(";")
                result[line[0]] = line[1]
    return result


def read_file_app(path):
    result = {}
    with open(path) as f:
        for l in f:
            if l.strip():
                line = l.strip().split(";")
                result[line[0].split("_")[0]] = line[1]
    return result


def count_positions(seq, aa="H"):
    positions = []
    for e, aa_seq in enumerate(seq):
        if len(aa) == 1:
            if aa == aa_seq:
                positions.append(e)
        else:
            for a in aa:
                if a == aa_seq:
                    positions.append(f"{a}_{e}")
    return positions


def save_result(path, organisms_species, organisms, sequences, organisms_hkr, organisms_ed, z_mean_appa,
                z_mean_total_appa, z_mean_appa_16plus,
                z_mean_appb, z_mean_total_appb, z_mean_appb_16plus, ncbi, tax_important, simply):
    with open(path, "w") as f:
        for species in organisms_species:
            lineage = ncbi.get_lineage(species)
            text_line = ""
            for lin in lineage:
                print(ncbi.get_rank([lin])[lin], tax_important)
                if "species" == ncbi.get_rank([lin])[lin]:
                    species2 = lin
                if simply:
                    if ncbi.get_rank([lin])[lin] in tax_important:
                        text_line = add_line_appa_appb(lin, ncbi, organisms, organisms_ed, organisms_hkr, text_line,
                                                       z_mean_appa, z_mean_appa_16plus, z_mean_appb, z_mean_appb_16plus,
                                                       z_mean_total_appa, z_mean_total_appb)
                else:
                    text_line = add_line_appa_appb(lin, ncbi, organisms, organisms_ed, organisms_hkr, text_line,
                                                   z_mean_appa, z_mean_appa_16plus, z_mean_appb, z_mean_appb_16plus,
                                                   z_mean_total_appa, z_mean_total_appb)
            if species2 is not None:
                text_line += str(sequences[species]) + "\t"
            species2 = None
            f.write(text_line + "\n")


def add_line_appa_appb(lin, ncbi, organisms, organisms_ed, organisms_hkr, text_line, z_mean_appa, z_mean_appa_16plus,
                       z_mean_appb, z_mean_appb_16plus, z_mean_total_appa, z_mean_total_appb):
    text_line += str(
        lin) + f"({ncbi.get_rank([lin])[lin]}, {ncbi.get_taxid_translator([lin])[lin]})" + "\t"
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms[lin].items()]) + "\t")
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms_hkr[lin].items()]) + "\t")
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms_ed[lin].items()]) + "\t")
    if len(z_mean_appa[lin]) == 0:
        text_line += "0"
    else:
        text_line += str(sum(z_mean_appa[lin]) / len(z_mean_appa[lin])) + "\t"
    if len(z_mean_total_appa[lin]) == 0:
        text_line += "0"
    else:
        text_line += str(sum(z_mean_total_appa[lin]) / len(z_mean_total_appa[lin])) + "\t"
    if len(z_mean_appa_16plus[lin]) == 0:
        text_line += "0"
    else:
        text_line += str(sum(z_mean_appa_16plus[lin]) / len(z_mean_appa_16plus[lin])) + "\t"
    if len(z_mean_appb[lin]) == 0:
        text_line += "0"
    else:
        text_line += str(sum(z_mean_appb[lin]) / len(z_mean_appb[lin])) + "\t"
    if len(z_mean_total_appb[lin]) == 0:
        text_line += "0"
    else:
        text_line += str(sum(z_mean_total_appb[lin]) / len(z_mean_total_appb[lin])) + "\t"
    if len(z_mean_appb_16plus[lin]) == 0:
        text_line += "0"
    else:
        text_line += str(sum(z_mean_appb_16plus[lin]) / len(z_mean_appb_16plus[lin])) + "\t"
    return text_line


def create_files(final_file, main_organisms, result_folder):
    tax_important = [
        'order',
        'family',
        'genus',
        "species",
        "subspecies",
    ]

    data = read_file(final_file)
    ncbi = NCBITaxa()
    sequences = {}
    for main, list_o in main_organisms.items():
        organisms_species = set()
        organisms = {}
        organisms_hkr = {}
        organisms_ed = {}
        z_mean = {}
        z_mean_16plus = {}
        organisms_total = {}
        organisms_hkr_total = {}
        organisms_hkr_16plus = {}
        organisms_ed_16plus = {}
        organisms_ed_total = {}
        z_mean_total = {}
        main_organism_tax_ids = list_o.keys()
        for e, res in enumerate(data.items()):
            tax_id, sequence = res
            sequences[tax_id] = sequence
            print(e)
            lineage = ncbi.get_lineage(tax_id)
            is_this = False
            for i in main_organism_tax_ids:
                if i in lineage:
                    is_this = True
            if is_this:
                positions, positions_ed, positions_ed_16plus, positions_ed_total, positions_hkr, positions_hkr_16plus, positions_hkr_total, positions_total = count_all_positions(
                    sequence)

                organisms_species.add(tax_id)
                for lin in lineage:
                    if lin not in organisms:
                        organisms[lin] = {}
                        organisms_hkr[lin] = {}
                        organisms_ed[lin] = {}
                        organisms_ed_16plus[lin] = {}
                        organisms_hkr_16plus[lin] = {}
                        z_mean[lin] = []
                        z_mean_16plus[lin] = []
                        organisms_total[lin] = {}
                        organisms_hkr_total[lin] = {}
                        organisms_ed_total[lin] = {}
                        z_mean_total[lin] = []
                    count_aa_positions(lin, organisms, organisms_ed, organisms_ed_16plus, organisms_ed_total,
                                       organisms_hkr,
                                       organisms_hkr_16plus, organisms_hkr_total, organisms_total, positions,
                                       positions_ed,
                                       positions_ed_16plus, positions_ed_total, positions_hkr, positions_hkr_16plus,
                                       positions_hkr_total, positions_total)

                    z_mean[lin].append(len(positions_hkr) - len(positions_ed))
                    z_mean_total[lin].append(len(positions_hkr_total) - len(positions_ed_total))
                    z_mean_16plus[lin].append(len(positions_hkr_16plus) - len(positions_ed_16plus))

        print(organisms)
        with open(os.path.join(result_folder, f"{main}_simpl.csv"), "w") as f:
            for species in organisms_species:
                lineage = ncbi.get_lineage(species)
                text_line = ""
                for lin in lineage:
                    print(ncbi.get_rank([lin])[lin], tax_important)
                    if ncbi.get_rank([lin])[lin] in tax_important:
                        text_line = add_line(lin, ncbi, organisms, organisms_ed, organisms_hkr, organisms_hkr_total,
                                             text_line, z_mean, z_mean_16plus, z_mean_total)
                text_line += str(sequences[species]) + "\t" + species
                f.write(text_line + "\n")
        with open(os.path.join(result_folder, f"{main}.csv"), "w") as f:
            for species in organisms_species:
                lineage = ncbi.get_lineage(species)
                text_line = ""
                for lin in lineage:
                    text_line = add_line(lin, ncbi, organisms, organisms_ed, organisms_hkr, organisms_hkr_total,
                                         text_line, z_mean, z_mean_16plus, z_mean_total)
                text_line += str(sequences[species]) + "\t" + species
                f.write(text_line + "\n")

    organisms_species = set()
    organisms = {}
    organisms_hkr = {}
    organisms_ed = {}
    z_mean_appa = {}
    z_mean_appb = {}

    organisms_hkr_16plus = {}
    organisms_ed_16plus = {}
    z_mean_appa_16plus = {}
    z_mean_appb_16plus = {}

    organisms_total = {}
    organisms_hkr_total = {}
    organisms_ed_total = {}
    z_mean_total_appa = {}
    z_mean_total_appb = {}

    sequences = {}
    for i in ["appa", "appb"]:
        data = read_file_app(os.path.join(result_folder, f"final_file_{i}_aligned.csv"))
        for e, res in enumerate(data.items()):
            tax_id, sequence = res
            if tax_id in sequences:
                sequences[tax_id] += ";" + i + ":" + sequence
            else:
                sequences[tax_id] = i + ":" + sequence
            if "canonical" not in tax_id:
                lineage = ncbi.get_lineage(tax_id)
                positions, positions_ed, positions_ed_16plus, positions_ed_total, positions_hkr, positions_hkr_16plus, positions_hkr_total, positions_total = count_all_positions(
                    sequence)

                organisms_species.add(tax_id)
                for lin in lineage:
                    if lin not in organisms:
                        organisms[lin] = {}
                        organisms_hkr[lin] = {}
                        organisms_ed[lin] = {}
                        organisms_total[lin] = {}
                        organisms_hkr_total[lin] = {}
                        organisms_ed_total[lin] = {}
                        organisms_hkr_16plus[lin] = {}
                        organisms_ed_16plus[lin] = {}

                        z_mean_appa_16plus[lin] = []
                        z_mean_appb_16plus[lin] = []
                        z_mean_appa[lin] = []
                        z_mean_appb[lin] = []
                        z_mean_total_appa[lin] = []
                        z_mean_total_appb[lin] = []
                    count_aa_positions(lin, organisms, organisms_ed, organisms_ed_16plus, organisms_ed_total,
                                       organisms_hkr, organisms_hkr_16plus, organisms_hkr_total, organisms_total,
                                       positions + [i], positions_ed + [i],
                                       positions_ed_16plus + [i],
                                       positions_ed_total + [i],
                                       positions_hkr + [i],
                                       positions_hkr_16plus + [i],
                                       positions_hkr_total + [i],
                                       positions_total + [i])

                    if i == "appa":
                        z_mean_appa[lin].append(len(positions_hkr) - len(positions_ed))
                        z_mean_total_appa[lin].append(len(positions_hkr_total) - len(positions_ed_total))
                        z_mean_appa_16plus[lin].append(len(positions_hkr_16plus) - len(positions_ed_16plus))
                    elif i == "appb":
                        z_mean_appb[lin].append(len(positions_hkr) - len(positions_ed))
                        z_mean_total_appb[lin].append(len(positions_hkr_total) - len(positions_ed_total))
                        z_mean_appb_16plus[lin].append(len(positions_hkr_16plus) - len(positions_ed_16plus))

    save_result(os.path.join(result_folder, f"actino_appa_appb.csv"), organisms_species, organisms, sequences,
                organisms_hkr, organisms_ed, z_mean_appa,
                z_mean_total_appa, z_mean_appa_16plus,
                z_mean_appb, z_mean_total_appb, z_mean_appb_16plus, ncbi, tax_important, simply=False)
    save_result(os.path.join(result_folder, f"actino_appa_appb_simpl.csv"), organisms_species, organisms, sequences,
                organisms_hkr, organisms_ed, z_mean_appa,
                z_mean_total_appa, z_mean_appa_16plus,
                z_mean_appb, z_mean_total_appb, z_mean_appb_16plus, ncbi, tax_important, simply=True)


def count_aa_positions(lin, organisms, organisms_ed, organisms_ed_16plus, organisms_ed_total, organisms_hkr,
                       organisms_hkr_16plus, organisms_hkr_total, organisms_total, positions, positions_ed,
                       positions_ed_16plus, positions_ed_total, positions_hkr, positions_hkr_16plus,
                       positions_hkr_total, positions_total):
    if tuple(positions) not in organisms[lin]:
        organisms[lin][tuple(positions)] = 0
    organisms[lin][tuple(positions)] += 1
    if tuple(positions_hkr) not in organisms_hkr[lin]:
        organisms_hkr[lin][tuple(positions_hkr)] = 0
    organisms_hkr[lin][tuple(positions_hkr)] += 1
    if tuple(positions_hkr_16plus) not in organisms_hkr_16plus[lin]:
        organisms_hkr_16plus[lin][tuple(positions_hkr_16plus)] = 0
    organisms_hkr_16plus[lin][tuple(positions_hkr_16plus)] += 1
    if tuple(positions_ed_16plus) not in organisms_ed_16plus[lin]:
        organisms_ed_16plus[lin][tuple(positions_ed_16plus)] = 0
    organisms_ed_16plus[lin][tuple(positions_ed_16plus)] += 1
    if tuple(positions_total) not in organisms_total[lin]:
        organisms_total[lin][tuple(positions_total)] = 0
    organisms_total[lin][tuple(positions_total)] += 1
    if tuple(positions_hkr_total) not in organisms_hkr_total[lin]:
        organisms_hkr_total[lin][tuple(positions_hkr_total)] = 0
    organisms_hkr_total[lin][tuple(positions_hkr_total)] += 1
    if tuple(positions_ed_total) not in organisms_ed_total[lin]:
        organisms_ed_total[lin][tuple(positions_ed_total)] = 0
    organisms_ed_total[lin][tuple(positions_ed_total)] += 1
    if tuple(positions_ed) not in organisms_ed[lin]:
        organisms_ed[lin][tuple(positions_ed)] = 0
    organisms_ed[lin][tuple(positions_ed)] += 1


def count_all_positions(sequence):
    positions = count_positions(sequence[:16])
    positions_hkr = count_positions(sequence[:16], aa="HKR")
    positions_ed = count_positions(sequence[:16], aa="ED")
    positions_hkr_16plus = count_positions(sequence[16:], aa="HKR")
    positions_ed_16plus = count_positions(sequence[16:], aa="ED")
    positions_total = count_positions(sequence)
    positions_hkr_total = count_positions(sequence, aa="HKR")
    positions_ed_total = count_positions(sequence, aa="ED")
    return positions, positions_ed, positions_ed_16plus, positions_ed_total, positions_hkr, positions_hkr_16plus, positions_hkr_total, positions_total


def add_line(lin, ncbi, organisms, organisms_ed, organisms_hkr, organisms_hkr_total, text_line, z_mean, z_mean_16plus,
             z_mean_total):
    text_line += str(
        lin) + f"({ncbi.get_rank([lin])[lin]}, {ncbi.get_taxid_translator([lin])[lin]})" + "\t"
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms[lin].items()]) + "\t")
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms_hkr[lin].items()]) + "\t")
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms_ed[lin].items()]) + "\t")
    text_line += str(';'.join([f"{i}:{j}" for i, j in organisms_hkr_total[lin].items()]) + "\t")
    text_line += str(sum(z_mean[lin]) / len(z_mean[lin])) + "\t"
    text_line += str(sum(z_mean_total[lin]) / len(z_mean_total[lin])) + "\t"
    text_line += str(sum(z_mean_16plus[lin]) / len(z_mean_16plus[lin])) + "\t"
    return text_line


@click.option('--final_file', default="./result/final_file_aligned.csv",
              help='Path to finale file from script "select_nonredundant_proteins.py".')
@click.option('--result_folder', default="./result",
              help='Path to folder with result files.".')
@click.option('--organisms', default={"7777": {7777: 'Chondrichthyes'},
                                      "7898": {7898: 'Actinopterygii', 7878: "Dipnomorpha",
                                               118072: "Coelacanthimorpha"},
                                      "40674": {40674: 'Mammalia'},
                                      "8782": {8782: 'Aves'},
                                      "1294634": {1294634: 'Crocodylia', 8459: 'Testudines', 8504: 'Lepidosauria'},
                                      '8292': {8292: 'Amphibia'}
                                      },
              help='''Dictionary in format:
              {'TaxIDGeneral:{TaxId:'OrganismName'}'}
              where based on TaxIDGeneral, there is created tree topology, organisms names will be in leaves of result tree.  
              if empty:
              '{"7777": {7777: 'Chondrichthyes'},
                 "7898": {7898: 'Actinopterygii', 7878: "Dipnomorpha", 118072: "Coelacanthimorpha"},
                 "40674": {40674: 'Mammalia'},
                 "8782": {8782: 'Aves'},
                 "1294634": {1294634: 'Crocodylia', 8459: 'Testudines', 8504: 'Lepidosauria'},
                 '8292': {8292: 'Amphibia'}
                }'
              ''')
def run(final_file, result_folder, organisms):
    if organisms:
        main_organisms = eval(organisms)
    else:
        main_organisms = \
            {"7777": {7777: 'Chondrichthyes'},
             "7898": {7898: 'Actinopterygii', 7878: "Dipnomorpha", 118072: "Coelacanthimorpha"},
             "40674": {40674: 'Mammalia'},
             "8782": {8782: 'Aves'},
             "1294634": {1294634: 'Crocodylia', 8459: 'Testudines', 8504: 'Lepidosauria'},
             '8292': {8292: 'Amphibia'}
             }
    create_files(final_file, main_organisms, result_folder)


if __name__ == "__main__":
    run()
