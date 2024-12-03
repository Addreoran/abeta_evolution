import PyQt5
from ete3 import NCBITaxa, TreeStyle, NodeStyle, AttrFace
from ete3.treeview import faces

from protein import Protein

def get_redundancy(proteins):
    new_proteins = []
    for protein in proteins:
        added = False
        for protein2 in proteins:
            if protein2 not in new_proteins:
                if protein != protein2:
                    if protein.seq == protein2.seq and protein.organism_id == protein2.organism_id:
                        print(protein.uniprot_id, protein2.uniprot_id)
                        protein.sequence_redundancy.append(protein2)
                        if protein not in new_proteins:
                            new_proteins.append(protein)
                        added = True
        if not added:
            new_proteins.append(protein)
    return new_proteins


def read_datas(file):
    file = open(file, "r")
    proteins = []
    for i in file.readlines():
        proteins.append(Protein(uniprot_id=i.split(";")[-1], organism_id=i.split(";")[0], seq=i.split(";")[1]))
    return proteins


def read_uniq_datas(file, positions):
    file = open(file, "r")
    proteins = {}
    empty = []
    for i in file.readlines():
        protein = Protein(i.split(";")[0], i.split(";")[2], i.split(";")[1], i.split(";")[-1])
        if i.split(";")[-1].strip() != "":
            if (i.split(";")[1], i.split(";")[-1]) in proteins.keys():
                if not proteins[(i.split(";")[1], i.split(";")[-1])].get_right_place_histidine(positions):
                    if protein.get_right_place_histidine(positions):
                        proteins[(i.split(";")[1], i.split(";")[-1])] = protein
            else:
                proteins[(i.split(";")[1], i.split(";")[-1])] = protein
        elif i.split(";")[-1].strip() == "":
            empty.append(protein)
        else:
            proteins[(i.split(";")[1], i.split(";")[-1])] = protein
    return list(proteins.values()) + empty

    # https://www.uniprot.org/uniprot/?query=A0A0S7JE63&columns=id,organism-id&format=tab


def get_class(organisms, ncbi):
    kingdoms = []
    classes = [40674, 8292, 8782, 7777, 7898, 118072, 8509, 1294634, 8459]

    # for i in organisms:
    #     ranks = ncbi.get_rank(ncbi.get_lineage(i))

    for org in organisms:
        is_king = False
        ranks = ncbi.get_lineage(org)
        # print(ranks)
        for cla in classes:
            if cla in ranks:
                if not (8509 in ranks or 8492 in ranks) and 32561 in ranks:
                    print("", org, ranks)
                    print(cla)

                kingdoms.append(cla)
                is_king = True
        if not is_king:
            print("brak", org, ncbi.get_rank(ncbi.get_lineage(org)))
        # if not (40674 in ranks or 8782 in ranks) and 32524 in ranks:
        #     print("owodniowiec", org, ranks)
    kingdoms = set(kingdoms)
    return kingdoms


def run(positions, type="general"):
    ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()
    # positions = {"H": 5}
    proteins = read_datas("../data/final_results/final_file.csv")
    # check_proteins("./dane_without_gene_redundancy.csv")

    # print(len(proteins))
    # proteins = get_redundancy(proteins)
    # proteins = read_uniq_datas("dane.csv", positions)
    # print(len(proteins))
    # proteins = [i for i in proteins if 40674 in ncbi.get_lineage(i.get_organism(ncbi=ncbi))]
    organisms = [i.get_organism(ncbi=ncbi, typ="diff") for i in proteins]
    print("select organisms", organisms)
    kingdoms = get_class(organisms, ncbi)
    print("select organisms", kingdoms)
    # organisms = set([i.get_organism(ncbi=ncbi, typ="general", selected=kingdoms) for i in proteins])
    # print("select organism tree", len(organisms))
    groups = {i: [] for i in kingdoms}
    groupped = 0
    nr = 15
    for i in [int(i) for i in set(organisms)]:
        lineage = ncbi.get_lineage(i)
        for group in groups.keys():
            if group in lineage:
                groups[group].append(i)

    organism_tree = ncbi.get_topology([int(i) for i in kingdoms])

    t = organism_tree
    print("plot tree")
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -45  # 0 degrees = 3 o'clock
    ts.arc_span = 90

    z_H = {}
    encount = {}
    print("dlugosci ", len(proteins))
    for i in proteins:
        lineage = ncbi.get_lineage(i.get_organism(ncbi=ncbi, typ="diff"))
        for line in lineage:
            if line not in z_H.keys():
                z_H[line] = {}
            for aa, pos in i.find_by_positions_histidine(positions).items():
                if aa not in z_H[line].keys():
                    z_H[line][aa] = pos
                else:
                    z_H[line][aa] += pos

            # if i.get_right_place_histidine(positions):
            #     z_H[i.get_organism(ncbi=ncbi)] += 1
            if line not in encount.keys():
                encount[line] = 1
            else:
                encount[line] += 1

    def layout(node):
        if node.is_leaf():
            # Add node name to laef nodes
            # N = AttrFace("name", fsize=100, fgcolor="black") #new
            N = AttrFace("name", fsize=60, fgcolor="black")

            faces.add_face_to_node(N, node, 0)
        if "weight" in node.features:
            # C = CircleFace(radius=node.weight, color="red", style="sphere")
            # # Let's make the sphere transparent
            # C.opacity = 0.3
            # # And place as a float face over the tree
            # faces.add_face_to_node(C, node, 0, position="float")
            pass

    nstyle = NodeStyle()
    # nstyle["hz_line_type"] = 1
    nstyle["hz_line_width"] = 20
    nstyle["hz_line_color"] = "blue"
    nstyle["vt_line_width"] = 20
    nstyle["vt_line_color"] = "blue"
    for n in t.traverse():
        suma_H = 0
        suma_all = 0

        # waga = 100 * suma_H / encount[n.taxid]
        # n.add_features(weight=waga, size=100)
        if len(n.get_leaf_names()) > 1:
            position_str = {}
            print(n.get_leaf_names())
            leaves = [int(i) for i in n.get_leaf_names()]
            new_pos = []
            # for z, j in suma_H_node.items():
            #     if z in leaves:
            #         for pos, g in j.items():
            #             if pos not in position_str:
            #                 position_str[pos] = g
            #             else:
            #                 position_str[pos] += g
            mapper = {"H06": "H6", "H213": "H13", "H314": "H14", "Q06": "Q6", "Q213": "Q13", "Q314": "Q14", }

            for z, j in z_H[n.taxid].items():
                new_pos.append(f"{mapper[z]}:{j}")
            position_str = "\n".join(new_pos)

            longNameFace = faces.TextFace(n.sci_name + "\n" + "All:" + str(encount[n.taxid]) + "\n" + str(position_str),
                                          fsize=60)

            n.add_face(longNameFace, column=0)
            n.set_style(nstyle)
    nstyle2 = NodeStyle()
    # nstyle2["hz_line_type"] = 1
    nstyle2["hz_line_width"] = 20
    nstyle2["hz_line_color"] = "red"
    nstyle2["vt_line_width"] = 20
    nstyle2["vt_line_color"] = "red"
    nstyle3 = NodeStyle()
    # nstyle3["hz_line_type"] = 1
    nstyle3["hz_line_width"] = 20
    nstyle3["hz_line_color"] = "violet"
    nstyle3["vt_line_width"] = 20
    nstyle3["vt_line_color"] = "violet"
    for i in t.get_leaves():
        tax_id = ncbi.get_taxid_translator([i.get_leaf_names()[0]])
        print(tax_id)
        position_str = []
        mapper = {"H06": "H6", "H213": "H13", "H314": "H14", "Q06": "Q6", "Q213": "Q13", "Q314": "Q14", }
        for z, j in z_H[i.taxid].items():
            position_str.append(f"{mapper[z]}:{j}")
        position_str = "\n".join(position_str)
        i.name = tax_id[int(i.get_leaf_names()[0])] + "\nAll:" + str(
            encount[int(i.get_leaf_names()[0])]) + "\n" + position_str
        if z_H[i.taxid]["H06"] == 0:
            i.set_style(nstyle2)
        elif z_H[i.taxid]["H06"] < encount[i.taxid] / 2:
            i.set_style(nstyle3)
        else:
            i.set_style(nstyle)
        i.dist = 5

    ts.layout_fn = layout
    ts.scale = 100  # new
    ts.show_leaf_name = False
    # ts.legend.add_face(TextFace("0.5 support"), column=1)
    t.render(f"../data/final_results/drzewo_blue.jpg", tree_style=ts)
    t.show(tree_style=ts)


def check_proteins(file):
    file = open(file, "r")
    proteins = []
    for i in file.readlines():
        print(i)
        proteins.append(i.split(";")[0])
    print("check redundancy", len(set(proteins)), len(proteins))


run({"H0": 5, "H2": 12, "H3": 13})
