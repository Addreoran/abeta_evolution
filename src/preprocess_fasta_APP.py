import os
import re
from multiprocessing import Pool

import requests

excluded_proteins = {"F6ZXB5", 'A0A8C9WM95', 'A0A6Q2YHZ8', 'A0A2I3MME2', 'A0A2I3N8B6', 'A0A3Q3GI80', 'C0JV68'}


def get_APP_protein(file_uniprot, file_app):
    file_app_new = open(file_app, "w")
    with open(file_uniprot) as f:
        is_APP = False
        for line in f.readlines():
            if line.strip():
                if line.startswith(">"):
                    is_app_re = re.search("\sGN=[Aa][pP][pP]\s", line)
                    if is_app_re:
                        is_APP = True
                    else:
                        print(is_APP, line)
                        is_APP = False
                if is_APP:
                    file_app_new.write(line)
    file_app_new.close()


def run_mafft(file, out):
    os.system(f'"/usr/bin/mafft"  --anysymbol --auto --clustalout --reorder "{file}" > "{out}"')


def search_APP_localisation(file_aln, file_out_aln, file_out_aln_excluded):
    sequence = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
    acc_human = "P05067"
    sequences = {}
    with open(file_aln) as f:
        # 1 otworzyć i spisać sekwencje
        for line in f.readlines():
            if line.strip() and "CLUSTAL" not in line:
                acc, sequence_acc = line.strip().split()
                if acc not in sequences:
                    sequences[acc] = sequence_acc
                else:
                    sequences[acc] += sequence_acc
                if acc_human in acc:
                    acc_human = acc
    # 2) znaleźć lokalizację
    pattern = r""
    for aa in sequence:
        pattern += fr"{aa}[-]*"
    res = re.search(pattern, sequences[acc_human])
    begin = res.start()
    end = res.end()
    sequences = {i: j[begin:end] for i, j in sequences.items()}
    # 3) zapisać tylko niezbędne alignmenty
    sequence_motif = "HDSGYEVHH"
    pattern2 = r""
    for aa in sequence_motif:
        pattern2 += fr"{aa}[-]*"
    res = re.search(pattern2, sequences[acc_human])
    begin = res.start()
    end = res.end()
    sequences_included = {i: j for i, j in sequences.items() if len(j[begin:end].replace("-", "")) > 5}
    sequences_excluded = {i: j for i, j in sequences.items() if len(j[begin:end].replace("-", "")) < 5}

    with open(file_out_aln, "w") as f:
        for acc, sequence_AB in sequences_included.items():
            f.write(f"{acc}\t{sequence_AB}\n")
    with open(file_out_aln_excluded, "w") as f:
        for acc, sequence_AB in sequences_excluded.items():
            f.write(f"{acc}\t{sequence_AB}\n")
    return sequences


def aln_to_fasta(file_fasta_all, file_fasta, sequences):
    headers = set()
    with open(file_fasta_all) as f:
        for line in f.readlines():
            if line.startswith(">"):
                headers.add(line)
            line = line.strip().split("\t")
    with open(file_fasta, "w") as f:
        for identificator, sequence in sequences.items():
            header_id = None
            for header in headers:
                if identificator in header:
                    if header_id is None:
                        header_id = header
                        break
            f.write(header_id)
            f.write(sequence.replace("-", ""))
            f.write("\n")


def get_excluded_names(file_fasta, file_excluded, file_aln, analise_file):
    acc_headers = {}
    organism_acc = {}
    with open(file_fasta) as f:
        for line in f.readlines():
            if line.startswith(">"):
                acc = line.split("|")[1]
                if acc not in acc_headers:
                    # >tr|A0A0A0MPX8|A0A0A0MPX8_FELCA Amyloid-beta A4 protein OS=Felis catus OX=9685 GN=APP PE=3 SV=2
                    organism = line.split("OX=")[-1].split("OS=")[-1].split()[0]
                    acc_headers[acc] = {"organism": organism,
                                        "taxonomy": line.split("OX=")[-1].split()[0],
                                        "name": line.split(" ", 1)[1].split("OS=")[0]}
                if organism not in organism_acc:
                    organism_acc[organism] = set()
                organism_acc[organism].add(acc)
    acc_excluded = set()
    with open(file_excluded) as f:
        for line in f.readlines():
            if line.strip():
                acc_excluded.add(line.split("|")[1])
    acc_aligned = set()
    with open(file_aln) as f:
        for line in f.readlines():
            if line.strip():
                acc_aligned.add(line.split("|")[1])
    with open(analise_file, "w") as f:
        for acc in acc_excluded:
            organism = acc_headers[acc]["organism"]
            organism_accs_icluded = organism_acc[organism].intersection(acc_aligned)
            print(acc_headers[acc])
            f.write(f"{acc}; {acc_headers[acc]['name']}; {organism}; {organism_accs_icluded}\n")


# import os
#
# # get the working directory
# base_path = os.getcwd()
#
# # this is the name of your file (with the 'py' extension)
# your_file_name = "imieninywybor.py"
#
# # Concatenate the file to your path
# full_path = os.path.join(base_path, your_file_name)
def get_uniprot(file_uniprot):
    import urllib.request
    urllib.request.urlretrieve(
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
        "../data/uniprot_trembl.fasta.gz")
    # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
    urllib.request.urlretrieve(
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
        "../data/uniprot_sprot.fasta.gz")
    # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    os.system(f"zcat ../data/uniprot_trembl.fasta.gz > {file_uniprot}")
    os.system(f"zcat ../data/uniprot_sprot.fasta.gz >> {file_uniprot}")
    # połączyć je
    pass


def get_fasta_of_aln(file_aln, fasta_all, fasta_ab):
    acc_headers = {}
    with open(fasta_all) as f:
        for line in f.readlines():
            if line.startswith(">"):
                acc = line.split("|")[1]
                if acc not in acc_headers:
                    # >tr|A0A0A0MPX8|A0A0A0MPX8_FELCA Amyloid-beta A4 protein OS=Felis catus OX=9685 GN=APP PE=3 SV=2
                    acc_headers[acc] = line
    seq_data = {}
    with open(file_aln) as f:
        for line in f.readlines():
            if line.strip():
                line = line.split()
                acc = line[0].split("|")[1]
                sequence = line[1].replace("-", "")
                seq_data[acc] = sequence
    if not os.path.exists("../data/fasta_ab/"):
        os.makedirs("../data/fasta_ab/")
    for acc, seq in seq_data.items():
        with open(f"../data/fasta_ab/{acc}.fasta", "w") as f:
            f.write(f">{acc_headers[acc]}\n")
            f.write(f"{seq}\n")


def run_jackhmmers(query_folder, result_folder, db_file):
    files = os.listdir(query_folder)
    if not os.path.exists(f"../data/{result_folder}"):
        os.makedirs(f"../data/{result_folder}")
    dataset = [(i.split(".")[0], db_file, result_folder) for i in files if "index" not in i]
    with Pool(200) as p:
        p.map(run_jackhmmer, dataset)


def run_jackhmmer(data):
    query_file, db_file, result_folder = data
    print(data)
    os.system(f"jackhmmer -A ../data/{result_folder}/jackhmmer_{query_file}.aln "
              f"-o ../data/{result_folder}/jackhmmer_{query_file}.txt --cpu 2 ../data/encode_seq/{query_file}.fasta {db_file}")
    pass


def read_fasta_accs(fasta_file):
    result = set()
    with open(fasta_file) as f:
        for line in f.readlines():
            if line.startswith(">"):
                if "|" in line:
                    result.add(line.split("|")[1])
                else:
                    result.add(line.strip().replace(">", ""))
    return result


def join_jackhmmer(jackhmmer_folder, output_fasta_file):
    files = os.listdir(jackhmmer_folder)
    files = [i for i in files if "aln" in i]
    new_acc = None
    new_similar = set()
    similar_graph = {}
    all_uniprots = set()
    fasta_file = open(output_fasta_file, "w")
    sequences = set()
    was = set()
    for file in files:
        with open(jackhmmer_folder + file) as f:
            for line in f.readlines():
                if line.strip():
                    if "STOCKHOLM" in line:
                        if new_acc is not None:
                            similar_graph[new_acc] = new_similar
                        new_acc = None
                        new_similar = set()
                    elif "GF ID" in line:
                        new_acc = line.split()[2].split("-")[0]
                    elif "#=GS" in line:
                        new_similar.add(line.split()[1].split("|")[1])
                    elif "#=GR" not in line and line.strip():
                        if "|" in line:
                            if line.split("|")[1] in new_similar:
                                line = line.split()
                                if line[0] not in was:
                                    was.add(line[0])
                                    fasta_file.write(f">{line[0]}\n")
                                    sequences.add(line[1].replace('-', ''))
                                    fasta_file.write(f"{line[1].replace('-', '')}\n")
                                    # all_uniprots.add(line.split()[1].split("|")[1])

    if new_acc is not None:
        similar_graph[new_acc] = new_similar
    fasta_file.close()
    print(len(all_uniprots), len(sequences))
    return similar_graph


def plot_pyvis(ab_swissprot, jackhmmer):
    from pyvis.network import Network
    pairs = {}
    for main_acc, minor_accs in jackhmmer.items():
        for minor_acc in minor_accs:
            if (main_acc, minor_acc) in pairs:
                pairs[(main_acc, minor_acc)] += 1
            elif (minor_acc, main_acc) in pairs:
                pairs[(minor_acc, main_acc)] += 1
            else:
                pairs[(main_acc, minor_acc)] = 1
    g = Network()
    added = set()
    print(len(pairs))
    for nodes, weight in pairs.items():
        color_1 = "blue"
        color_2 = "red"
        if nodes[0] not in added:
            if nodes[0] in ab_swissprot:
                g.add_node(n_id=nodes[0], label=nodes[0], color=color_2)
            else:
                g.add_node(n_id=nodes[0], label=nodes[0], color=color_1)
        if nodes[1] not in added:
            if nodes[1] in ab_swissprot:
                g.add_node(n_id=nodes[1], label=nodes[1], color=color_2)
            else:
                g.add_node(n_id=nodes[1], label=nodes[1], color=color_1)
        # g.add_node(1, color='red')
        added.add(nodes[0])
        added.add(nodes[1])
        g.add_edge(nodes[0], nodes[1],
                   # value=weight,
                   # title=str(weight)
                   )  # weight 42
    # g.save_graph("../data/nx.html")
    print(len(added))
    g.show_buttons(filter_=['physics'])

    g.show('nx.html', notebook=False)


def get_unique_sequences(file_aln, fasta_unique_ab):
    seq_uniq = {}
    with open(file_aln) as f:
        for line in f.readlines():
            if line.strip():
                line = line.strip().split()
                acc = line[0].split("|")[1]
                seq = line[1].replace("-", "")
            if seq not in seq_uniq:
                seq_uniq[seq] = set()
            seq_uniq[seq].add(acc)
    histidin_6 = 0
    with open(fasta_unique_ab, "w") as f:
        for seq, accs in seq_uniq.items():
            if len(accs.intersection(excluded_proteins)) == 0:
                f.write(f">{','.join(list(accs))}\n")
                f.write(seq + "\n")
                # print(seq, accs, accs.intersection(excluded_proteins))
                if seq[5] == "H":
                    histidin_6 += 1
                else:
                    print(seq, accs, accs.intersection(excluded_proteins))
    print(histidin_6)

    return seq_uniq


def encode_sequences(fasta_file, file_path, encode_fasta_file):
    if not os.path.exists(file_path):
        os.mkdir(file_path)
    ids_proteins = {}
    id = 0
    seq = None
    accs = None
    new = open(encode_fasta_file, "w")
    with open(fasta_file) as f:
        for line in f.readlines():
            # pass
            # with open(file_path, "w") as f2:
            if line.strip() and line.startswith(">"):
                if accs is not None:
                    with open(f"{file_path}/accs_{id}.fasta", "w") as f2:
                        f2.write(f">{id}\n")
                        f2.write(seq + "\n")
                        new.write(f">{id}\n")
                        new.write(seq + "\n")
                        ids_proteins[id] = accs
                        id += 1
                accs = line.replace(">", "")
            else:
                seq = line.strip()
    if accs is not None:
        with open(f"{file_path}/accs_{id}.fasta", "w") as f2:
            f2.write(f">{id}\n")
            f2.write(seq + "\n")
            new.write(f">{id}\n")
            new.write(seq + "\n")
            ids_proteins[id] = accs
            id += 1
    new.close()
    with open(f"{file_path}/index.txt", "w") as f:
        for ids, accs in ids_proteins.items():
            f.write(f"{ids}\t{accs}\n")


def encode_result(fasta_file, encode_fasta_file):
    ids_proteins = {}
    id = 0
    seq = None
    accs = None
    sequences = {}
    new = open(encode_fasta_file, "w")
    with open(fasta_file) as f:
        for line in f.readlines():
            if line.strip() and line.startswith(">"):
                if accs is not None:
                    if seq not in sequences:
                        new.write(f">result_{id}\n")
                        new.write(seq + "\n")
                        sequences[seq] = id
                        id += 1
                        ids_proteins[sequences[seq]] = set()
                    ids_proteins[sequences[seq]].add(accs)
                accs = line.strip().replace(">", "")
            else:
                seq = line.strip()
    if accs is not None:
        if seq not in sequences:
            new.write(f">result_{id}\n")
            new.write(seq + "\n")
            sequences[seq] = id
            id += 1
            ids_proteins[sequences[seq]] = set()
        ids_proteins[sequences[seq]].add(accs)
    new.close()
    with open(f"../data/index_result.txt", "w") as f:
        for ids, accs in ids_proteins.items():
            f.write(f"{ids}\t{','.join(list(accs))}\n")


def download_fasta(file_list):
    # if not os.path.exists("../data/organism/"):
    #     os.mkdir("../data/organism/")
    acc = set()
    for file in file_list:
        with open(file) as f:
            for line in f.readlines():
                if line.strip():
                    line = line.split("\t")
                    uniprot_list = line[1].split(",")
                    for uniprot in uniprot_list:
                        if "|" in uniprot:
                            acc.add(uniprot.strip().split("|")[1])
                        else:
                            acc.add(uniprot.strip())
    with open("../data/after_jackhmmer_total_sequences.fasta", "w") as f:
        for uniprot in acc:
            req = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta")
            f.write(req.text)


if __name__ == "__main__":
    # # 1) Pobranie białek powstających z genu APP
    # # https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28gene%3AAPP%29+AND+%28protein_name%3AAmyloid-beta%29%29
    # # 2) Usunięcie białek, które nie powstają z APP
    # get_APP_protein(file_uniprot="../data/uniprot_APP.fasta", file_app="../data/uniprot_only_APP.fasta")
    # # 3) szukamy fragmentów abeta
    # # 3A) mafft
    # run_mafft("../data/uniprot_only_APP.fasta", "../data/alignment_APP.aln")
    # # 3b) lokalizacja abety i zapisanie do pliku
    # sequences = search_APP_localisation(file_aln="../data/alignment_APP.aln", file_out_aln="../data/alignment_AB.aln",
    #                                     file_out_aln_excluded="../data/alignment_AB_excluded.aln")
    # # 3c) zapisanie fasta do pliku
    # aln_to_fasta(file_fasta_all="../data/uniprot_only_APP.fasta",
    #              file_fasta="../data/uniprot_AB.fasta",
    #              sequences=sequences)
    # # 3d) sprawdzić czy wykluczenia są ok
    # # - pobrać nazwy białek i orgainzmy,
    # # jeśli takie samo już jest w moim zestawie białek,
    # # to nie wyrzucam, jeśli nie ma, to sprawdzam czy na pewno powinno zostać wyrzucone
    # get_excluded_names(file_fasta="../data/uniprot_only_APP.fasta",
    #                    file_excluded="../data/alignment_AB_excluded.aln",
    #                    file_aln="../data/alignment_AB.aln",
    #                    analise_file="../data/excluded_acc_analyse.csv")
    # # 4) jackhmmer
    # # 4a) pobranie bazy trembl+sp
    # # get_uniprot("../data/uniprot.fasta")
    # # 4b) pobranie fasta z align
    # get_fasta_of_aln(file_aln="../data/alignment_AB.aln",
    #                  fasta_all="../data/uniprot_APP.fasta",
    #                  fasta_ab="../data/AB.fasta")
    # sequences = get_unique_sequences(file_aln="../data/alignment_AB.aln",
    #                                  fasta_unique_ab="../data/AB_uniq.fasta")
    # encode_sequences(file_path="../data/encode_seq/",
    #                  fasta_file="../data/AB_uniq.fasta",
    #                  encode_fasta_file="../data/encode_AB.fasta")
    # # 4c) puszczenie jackhmmer z input jako "../data/uniprot_AB.fasta"
    # # run_jackhmmers(query_folder="../data/encode_seq/",
    # #                result_folder="jackhmmer_encode",
    # #                db_file="../data/uniprot.fasta")
    # # # 4d) połączenie wyników jackhmmer
    # ad_fasta = read_fasta_accs("../data/encode_AB.fasta")
    # jackhmmer_result = join_jackhmmer("../data/jackhmmer_encode/",
    #                                   "../data/data_total_jackhmmer_sequences_encode.fasta")
    #
    # encode_result(fasta_file="../data/data_total_jackhmmer_sequences_encode.fasta",
    #               encode_fasta_file="../data/data_total_jackhmmer_sequences_encode_encode.fasta")
    # os.system("cat ../data/data_total_jackhmmer_sequences_encode_encode.fasta >../data/after_jackhmmer.fasta")
    # os.system("cat ../data/encode_AB.fasta >>../data/after_jackhmmer.fasta")
    # # plot_pyvis(ad_fasta, jackhmmer_result)
    # # # 5) mafft po raz drugi z usunięciem braków w alignmencie lub inne białka
    # run_mafft("../data/after_jackhmmer.fasta",
    #           "../data/after_jackhmmer.aln")
    # # znaleźć sekwencje białek i podzielić na organizmy
    download_fasta([
        "../data/encode_seq/index.txt",
        "../data/index_result.txt"
    ])
    run_mafft("../data/after_jackhmmer_total_sequences.fasta",
              "../data/after_jackhmmer_total_sequences.aln")
    sequences = search_APP_localisation(file_aln="../data/after_jackhmmer_total_sequences.aln",
                                        file_out_aln="../data/after_jackhmmer_total_sequences_AB.aln",
                                        file_out_aln_excluded="../data/after_jackhmmer_total_sequences_AB_excluded.aln")
    # 6) usunięcie redundancji
