import os

from ete3 import NCBITaxa, TreeStyle, AttrFace
import PyQt5

from src.protein import Protein


def read_H_position(path):
    result = {}
    with open(path) as f:
        header = f.readline().strip().split(";")
        for line in f:

            if line.strip() and "taxonomy" not in line:
                cut_line = line.strip().split(";")
                result[int(cut_line[0])] = {i: cut_line[e + 1] for e, i in enumerate(header[1:])}
    return result


def plot_tree(main_organisms, H_positions_file):
    global node
    H_positions = read_H_position(H_positions_file)
    main_organisms_ids = main_organisms.keys()
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(main_organisms_ids)
    for i in tree.traverse():
        children_node = i.children
        if children_node:
            if len(children_node[0].children) > len(children_node[1].children):
                i.swap_children()
    colors = {"All": "black",
              "many": "green",
              "small": "orange",
              "lack": "red"}
    selected_positions = {int(i): {"All": 0, "H06": 0, "H213": 0, "H314": 0} for i in main_organisms.keys()}
    for i, j in main_organisms.items():
        for tax_id in j.keys():
            selected_positions[int(i)]["All"] += int(H_positions[int(tax_id)].get("All", 0))
            selected_positions[int(i)]["H06"] += int(H_positions[int(tax_id)].get("H06", 0))
            selected_positions[int(i)]["H213"] += int(H_positions[int(tax_id)].get("H213", 0))
            selected_positions[int(i)]["H314"] += int(H_positions[int(tax_id)].get("H314", 0))
        pass
    for node in tree.iter_leaves():
        node.add_feature("All", "All: " + str(selected_positions[int(node.name)].get("All", "empty")))
        node.add_feature("H06", "H6: " + str(selected_positions[int(node.name)]["H06"]))
        node.add_feature("H13", "H13: " + str(selected_positions[int(node.name)]["H213"]))
        node.add_feature("H14", "H14: " + str(selected_positions[int(node.name)]["H314"]))
        all = int(selected_positions[int(node.name)].get("All", "0"))
        if int(selected_positions[int(node.name)]["H06"]) == all:
            H06_color = colors["many"]
        elif int(selected_positions[int(node.name)]["H06"]) >= all * 0.75:
            H06_color = colors["small"]
        elif int(selected_positions[int(node.name)]["H06"]) < all * 0.75:
            H06_color = colors["lack"]

        if int(selected_positions[int(node.name)]["H213"]) == all:
            H13_color = colors["many"]
        elif int(selected_positions[int(node.name)]["H213"]) >= all * 0.75:
            H13_color = colors["small"]
        elif int(selected_positions[int(node.name)]["H213"]) < all * 0.75:
            H13_color = colors["lack"]

        if int(selected_positions[int(node.name)]["H314"]) == all:
            H14_color = colors["many"]
        elif int(selected_positions[int(node.name)]["H314"]) >= all * 0.75:
            H14_color = colors["small"]
        elif int(selected_positions[int(node.name)]["H314"]) < all * 0.75:
            H14_color = colors["lack"]
        N = AttrFace("sci_name", fgcolor="black")
        N.bold = True

        all = AttrFace("All", fsize=8, fgcolor="black")
        H06 = AttrFace("H06", fsize=8, fgcolor=H06_color)
        H13 = AttrFace("H13", fsize=8, fgcolor=H13_color)
        H14 = AttrFace("H14", fsize=8, fgcolor=H14_color)
        node.sci_name = "\n"
        for names in main_organisms[str(node.name)].values():
            node.sci_name += names + "\n"
        node.sci_name = node.sci_name[:-1]
        N.margin_top = N.margin_bottom = N.margin_left = 4.0
        all.margin_top = all.margin_bottom = all.margin_left = 4.0
        node.add_face(N, 1, position='aligned')
        node.add_face(all, 1, position='aligned')
        node.add_face(H06, 1, position='aligned')
        node.add_face(H13, 1, position='aligned')
        node.add_face(H14, 1, position='aligned')
        node.name = " "
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.draw_guiding_lines = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    tree.show(tree_style=ts)


def read_data(file):
    file = open(file, "r")
    proteins = []
    for i in file.readlines():
        proteins.append(Protein(uniprot_id=i.split(";")[-1], organism_id=i.split(";")[0], seq=i.split(";")[-1]))
    return proteins


def get_class(organisms, ncbi, classes):
    kingdoms = []
    for org in organisms:
        is_king = False
        ranks = ncbi.get_lineage(org)
        for cla in classes:
            if cla in ranks:
                if not (8509 in ranks or 8492 in ranks) and 32561 in ranks:
                    print("", org, ranks)
                    print(cla)
                kingdoms.append(cla)
                is_king = True
        if not is_king:
            print("lack", org, ncbi.get_rank(ncbi.get_lineage(org)))
    kingdoms = set(kingdoms)
    return kingdoms


def get_H(final_file, H_positions_file, main_organisms, positions, ncbi):
    proteins = read_data(final_file)
    organisms = [i.get_organism(ncbi=ncbi, typ="diff") for i in proteins]
    classes = []
    for tax, group in main_organisms.items():
        classes += list(group.keys())
    kingdoms = get_class(organisms, ncbi, classes)
    groups = {i: [] for i in kingdoms}
    for i in [int(i) for i in set(organisms)]:
        lineage = ncbi.get_lineage(i)
        for group in groups.keys():
            if group in lineage:
                groups[group].append(i)
    z_H = {}
    encount = {}
    positions_values = set()
    for i in proteins:
        lineage = ncbi.get_lineage(i.get_organism(ncbi=ncbi, typ="diff"))
        for line in lineage:
            if line not in z_H.keys():
                z_H[line] = {}
            for aa, pos in i.find_by_positions_histidine(positions).items():
                if aa not in z_H[line].keys():
                    z_H[line][aa] = pos
                    positions_values.add(aa)
                else:
                    z_H[line][aa] += pos
            if line not in encount.keys():
                encount[line] = 1
            else:
                encount[line] += 1
    positions_values = list(positions_values)
    with open(H_positions_file, "a") as f:
        f.write("taxonomy_id;All;")
        f.write(';'.join(positions_values))
        f.write("\n")
        for line, data in z_H.items():
            f.write(f"{line};{encount[line]};")
            for pos in positions_values:
                f.write(f"{data[pos]};")
            f.write("\n")


import click


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
@click.option('--aa_positions', default={"H0": 5, "H2": 12, "H3": 13},
              help='''Dictionary in format:
              {'HPos:NumericPos'}
              if empty:
              '{"H0": 5, "H2": 12, "H3": 13}'
              ''')
@click.option('--final_file', default="./result/final_file_aligned.csv",
              help='Path to file with all TaxIds and peptides. File generated with script select_non_redundant_proteins.py')
@click.option('--H_positions_file', default="./result/H_positions_file.csv",
              help='Path to file that will contain positions of analised amino acids.')
def run(organisms=None, aa_positions=None, H_positions_file=None, final_file="./result/final_file_aligned.csv"):
    ncbi = NCBITaxa()

    if H_positions_file is None:
        H_positions_file = "./result/tmp_H_pos.csv"
    if not organisms:
        main_organisms = \
            {"7777": {7777: 'Chondrichthyes'},
             "7898": {7898: 'Actinopterygii', 7878: "Dipnomorpha", 118072: "Coelacanthimorpha"},
             "40674": {40674: 'Mammalia'},
             "8782": {8782: 'Aves'},
             "1294634": {1294634: 'Crocodylia', 8459: 'Testudines', 8504: 'Lepidosauria'},
             '8292': {8292: 'Amphibia'}
             }
    else:
        main_organisms = eval(organisms)
    if not aa_positions:
        aa_positions = {"H0": 5, "H2": 12, "H3": 13}
    else:
        aa_positions = eval(aa_positions)
    get_H(final_file, H_positions_file, main_organisms, aa_positions, ncbi)
    plot_tree(main_organisms, H_positions_file)
    os.remove(H_positions_file)


if __name__ == "__main__":
    run()
