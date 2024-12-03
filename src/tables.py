from ete3 import NCBITaxa

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


# def read_datas(file):
#     file = open(file, "r")
#     proteins = []
#     for i in file.readlines():
#         proteins.append(Protein(uniprot_id=i.split(";")[-1], organism_id=i.split(";")[0], seq=i.split(";")[1]))
#     return proteins

def prep_data(data, positions, groups):
    results = {}
    ncbi = NCBITaxa()
    # (group, aa, pos)={prot:org}
    for prot in data:
        for pos in positions.values():
            #             group=...
            group = prot.get_organism(ncbi, groups)
            if group in groups:
                # aa=...
                aa = prot.seq[pos]
                if aa != "H":
                    # print(ncbi.get_taxid_translator([str(group)]))
                    # print(ncbi.get_taxid_translator([prot.get_organism(ncbi, typ="diff")]))
                    # print(prot.get_organism(ncbi, typ="diff"))
                    # print(ncbi.get_taxid_translator(
                    #     [prot.get_organism(
                    #         ncbi,
                    #         typ="diff")])[
                    #           int(prot.get_organism(
                    #               ncbi,
                    #               typ="diff"))])
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
        f.write("tax_group;new amino acid;position;organism_no;proteins_info\n"
                "")
        for info, proteins in data.items():
            protein_data = [f"{i.uniprot_id}({j})" for i, j in set(proteins)]

            f.write(f"{info[0]};{info[1]};{info[2] + 1};"
                    f"{len(protein_data)}"
                    f";{','.join(protein_data)};"
                    f"\n")
    pass


groups = [7777, 7898, 118072, 40674, 8782, 8509, 1294634, 8459, 8292]
positions = {"H0": 5, "H2": 12, "H3": 13}
data = read_file("../data/final_results/final_file.csv")
res = prep_data(data, positions, groups)
save_table(res, "../data/final_results/prec_table.csv")
