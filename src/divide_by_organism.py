from ete3 import NCBITaxa

from protein import Protein


# todo:fix

def get_groups():
    groups = dict(
        fish=(7777, 7898, 118072),
        mammals=(40674),
        aves=(8782),
        reptiles=(8509, 1294634, 8459),
        amphibians=(8292)
    )
    return groups


def read_file(file):
    file = open(file, "r")
    proteins = []
    for i in file.readlines():
        print(i)
        protein = Protein(uniprot_id=i.strip().split(";")[-1], organism_id=i.strip().split(";")[0].replace('\n', ''),
                          seq=i.split(";")[1])
        protein.line = i
        proteins.append(protein)
    return proteins


def get_proteins_from_group(proteins, group):
    ncbi = NCBITaxa()
    proteins2 = [i for i in proteins if i.get_organism(ncbi, group) in group]
    return proteins2


def get_fish_with_Q(proteins, position, aa, taxId):
    ncbi = NCBITaxa()
    proteins_new = [i for i in proteins if i.get_organism(ncbi, typ="general", selected=taxId) in taxId]
    # print("org", [i.get_organism(ncbi, typ="general", selected=taxId) for i in proteins])
    # print(proteins_new)

    if aa:
        proteins_new = [i for i in proteins_new if i.seq[position] == aa]
    else:
        proteins_new = [i for i in proteins_new]
    return proteins_new


def save_group(proteins, file):
    with open(file, "w") as f:
        for p in proteins:
            f.write(p.line)


def run():
    proteins = read_file("../data/final_file.csv")
    groups = get_groups()
    proteins2 = get_fish_with_Q(proteins, 5, "Q", taxId=[7898])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_fish_Q.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[7777])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_Chondrichtlhyes.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[7898])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_Actinopterygii.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[8782])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_Birds.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[8292])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_Amphibians.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[8509, 1294634, 8459])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_reptiles.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[7777, 7898, 118072])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_fish.csv")
    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[7898])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_Coelacanthimorpha.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[40674])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_mammals.csv")

    proteins2 = get_fish_with_Q(proteins, 5, "", taxId=[7776])
    save_group(proteins2, f"../data/final_results/divided_groups_without_gene_redundancy_Gnathostomata.csv")

    # for name, group in groups.items():
    #     proteins2 = get_proteins_from_group(proteins, group)
    #     save_group(proteins2, f"divided_groups_without_gene_redundancy_{name}.csv")


run()
