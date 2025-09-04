from ete3 import NCBITaxa
import click
from protein import Protein


def read_file(file):
    file = open(file, "r")
    proteins = []
    for i in file.readlines():
        print("prot", i.split(";")[-1])
        protein = Protein(uniprot_id=i.strip().split(";")[-1], organism_id=i.strip().split(";")[0].replace('\n', ''),
                          seq=i.split(";")[1])
        protein.line = i
        proteins.append(protein)
    return proteins


def prep_data(data, positions, groups):
    results = {}
    ncbi = NCBITaxa()
    for prot in data:
        for pos in positions:
            group = prot.get_organism(ncbi, groups)
            if group in groups:
                aa = prot.seq[pos]
                if aa != "H":
                    if (ncbi.get_taxid_translator([str(group)])[group], aa, pos) not in results:
                        results[(ncbi.get_taxid_translator([str(group)])[group], aa, pos)] = [(prot,
                                                                                               ncbi.get_taxid_translator(
                                                                                                   [prot.get_organism(
                                                                                                       ncbi,
                                                                                                       typ="diff")])[
                                                                                                   int(prot.get_organism(
                                                                                                       ncbi,
                                                                                                       typ="diff"))])]

                    else:
                        results[(ncbi.get_taxid_translator([str(group)])[group], aa, pos)].append(
                            (prot, ncbi.get_taxid_translator([prot.get_organism(ncbi, typ="diff")])[int(
                                prot.get_organism(ncbi, typ="diff"))]))
    return results


def save_table(data, file):
    with open(file, "w") as f:
        f.write("tax_group;new amino acid;position;organism_no;proteins_info\n")
        for info, proteins in data.items():
            protein_data = [f"{i.uniprot_id}({j})" for i, j in set(proteins)]
            f.write(f"{info[0]};{info[1]};{info[2] + 1};"
                    f"{len(protein_data)}"
                    f";{','.join(protein_data)};"
                    f"\n")


@click.option('--final_file', default="./result/final_file_aligned.csv",
              help='Path to finale file from script "select_nonredundant_proteins.py".')
@click.option('--TaxIds',
              default=[40674, 8292, 8782, 7777, 7898, 118072, 1294634, 8459, 8504, 7878, 8509, 1294634, 8459, 32561],
              help='List of Taxonomy IDs to calculate statistics.')
@click.option('--aminoacids_positions', default=[5, 12, 13],
              help='Positions and names of aminoacids, that should be in table.')
@click.option('--result_file', default="./result/prec_table.csv",
              help='Path to save result table')

def run(final_file, TaxIds, aminoacids_positions, result_file):
    data = read_file(final_file)
    res = prep_data(data, aminoacids_positions, TaxIds)
    save_table(res, result_file)
