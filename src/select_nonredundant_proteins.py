import os
import re
import requests
from Levenshtein import distance
import click
from ete3 import NCBITaxa
import shutil
from bs4 import BeautifulSoup


def run_mafft(file, out):
    os.system(f'"/usr/bin/mafft"  --anysymbol --auto --clustalout --reorder "{file}" > "{out}"')


def search_APP_localisation(file_aln, file_out_aln, file_out_aln_excluded):
    sequence = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
    acc_human = "P05067"
    sequences = {}
    with open(file_aln) as f:
        for line in f.readlines():
            if line.strip() and "CLUSTAL" not in line:
                acc, sequence_acc = line.strip().split()
                if acc not in sequences:
                    sequences[acc] = sequence_acc
                else:
                    sequences[acc] += sequence_acc
                if acc_human in acc:
                    acc_human = acc
    pattern = r""
    for aa in sequence:
        pattern += fr"{aa}[-]*"
    res = re.search(pattern, sequences[acc_human])
    begin = res.start()
    end = res.end()
    sequences = {i: j[begin:end] for i, j in sequences.items()}
    sequence_motif = "HDSGYEVHH"
    pattern2 = r""
    for aa in sequence_motif:
        pattern2 += fr"{aa}[-]*"
    res = re.search(pattern2, sequences[acc_human])
    begin = res.start()
    end = res.end()
    sequences_included = {i: j for i, j in sequences.items() if len(j[begin:end].replace("-", "")) >= 5}
    sequences_excluded = {i: j for i, j in sequences.items() if len(j[begin:end].replace("-", "")) < 5}

    with open(file_out_aln, "w") as f:
        for acc, sequence_AB in sequences_included.items():
            f.write(f"{acc}\t{sequence_AB}\n")
    with open(file_out_aln_excluded, "w") as f:
        for acc, sequence_AB in sequences_excluded.items():
            f.write(f"{acc}\t{sequence_AB}\n")
    return sequences


def get_organisms_from_file(file):
    result = {}
    with open(file) as f:
        for line in f:
            if line.strip():
                acc, tax_id = line.strip().split(";")
                result[acc] = tax_id
    return result


def get_organism_by_refseq(refseq_id):
    ncbi = NCBITaxa()
    name2taxid = ncbi.get_name_translator([refseq_id])
    print(name2taxid)
    if len(name2taxid.get(refseq_id, {})) != 1:
        print(f"TaxId not found. Give tax id of {refseq_id}:")
        return input()
    return name2taxid[refseq_id][0]


def get_tax_uniprot(uniprot):
    utl = f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"
    req = requests.get(utl)
    for line in req.text.split("\n"):
        if line.startswith(">"):
            return line.split("OX=")[-1].split()[0]


def get_organisms_uniprot(file_aln_uniprot, fasta_file_uniprot, fasta_sequences=None, ox_sets=None):
    acc_ox = {}
    ox = None
    acc = None
    sequences = ""
    if fasta_sequences is None:
        fasta_sequences = {}
    if ox_sets is None:
        ox_sets = {}
    with open(fasta_file_uniprot) as f:
        for line in f.readlines():
            if line.startswith(">"):
                if acc is not None:
                    fasta_sequences[acc] = sequences
                    acc_ox[acc] = ox
                    ox = 0
                acc = line.split()[0].replace(">", "").split("|")[1]
                if "OX" in line:
                    tax = line.split("OX=")[-1].rsplit(" ", 2)[0]
                    if len(tax.split()) == 1 and not acc.startswith("UPI"):
                        ox = get_tax_uniprot(acc)
                    else:
                        ox = line.split("OX=")[-1].split()[0]
                sequences = line
            else:
                sequences += line
    if acc is not None:
        fasta_sequences[acc] = sequences
        acc_ox[acc] = ox
        ox = 0
    with open(file_aln_uniprot) as f2:
        for line in f2:
            if line.strip():
                acc = line.split()[0].replace(">", "").split("|")[1]
                if acc_ox.get(acc, "0") not in ox_sets:
                    ox_sets[acc_ox.get(acc, "0")] = set()
                ox_sets[acc_ox.get(acc, "0")].add(acc)
    return ox_sets, fasta_sequences


def get_organisms_refseq(file_aln_refseq, fasta_file_refseq, organism_file, fasta_sequences=None, ox_sets=None):
    if fasta_sequences is None:
        fasta_sequences = {}
    if ox_sets is None:
        ox_sets = {}
    if os.path.exists(organism_file):
        organism_refseq = get_organisms_from_file(organism_file)
    else:
        organism_refseq = {}
    acc = None
    ox = 0
    acc_ox = {}
    ox_tmp = {}
    sequences = ""
    with open(fasta_file_refseq) as f:
        for line in f.readlines():
            if line.startswith(">"):
                if acc is not None:
                    fasta_sequences[acc] = sequences
                    acc_ox[acc] = ox
                    ox = 0
                acc = line.split()[0].replace(">", "")
                if acc not in organism_refseq:
                    try:
                        if line.split("[")[-1].split("]")[0] in ox_tmp:
                            ox = ox_tmp[line.split("[")[-1].split("]")[0]]
                        else:
                            ox = get_organism_by_refseq(line.split("[")[-1].split("]")[0])
                            ox_tmp[line.split("[")[-1].split("]")[0]] = ox
                    except Exception as e:
                        print(f"{e}. TaxId not found. Give me TaxId of {acc}:")
                        ox = input().strip()
                    with open(organism_file, "a") as o:
                        o.write(f"{acc};{ox}\n")
                else:

                    ox = organism_refseq[acc]
                sequences = line
            else:
                sequences += line
    if acc is not None:
        fasta_sequences[acc] = sequences
        acc_ox[acc] = ox
        ox = 0
    acc_required = set()
    with open(file_aln_refseq) as f3:
        for line in f3:
            if line.strip():
                acc = line.split()[0].replace(">", "")
                acc_required.add(acc)
                if acc_ox.get(acc, "0") not in ox_sets:
                    ox_sets[acc_ox.get(acc, "0")] = set()
                ox_sets[acc_ox.get(acc, "0")].add(acc)
    return ox_sets, fasta_sequences


def divide_by_organisms(out_folder, ox_sets, fasta_sequences):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    for ox, accs in ox_sets.items():
        with open(f"{out_folder}/{ox}.fasta", "w") as f3:
            for acc in accs:
                f3.write(fasta_sequences[acc])


def update_organisms(fasta_sequences, ox_sets, out_folder, canonical_sequence):
    for ox, accs in ox_sets.items():
        with open(f"{out_folder}/{ox}.fasta", "w") as f:
            f.write(canonical_sequence)
            for acc in accs:
                if acc in fasta_sequences:
                    f.write(fasta_sequences[acc])


def make_mafft_per_organism(folder):
    files = [i for i in os.listdir(folder) if "fasta" in i]
    for file in files:
        if "index" not in file:
            run_mafft(folder + file, folder + file.split(".")[0] + ".aln")


def encode_mafft_find_amyloid_per_organism(folder, result_file, organims_with_isoforms):
    # {folder}index.csv
    set_id = 0
    sequence = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
    acc_human = "canonical_seq"
    pattern = r""
    index_accs = {}
    for aa in sequence:
        pattern += fr"{aa}[-]*"
    pattern = pattern[:-4]
    final_file = open(result_file, "w")
    with open(organims_with_isoforms, "w") as o:
        for file in [i for i in os.listdir(folder) if "aln" in i and "encoded" not in i]:
            with open(folder + file) as f:
                sequences = {}
                for line in f.readlines():
                    if line.strip() and "CLUSTAL" not in line and "*" not in line and len(line.strip().split()) == 2 and \
                            line.strip().split()[0] not in [".", "*", ":", "::"]:
                        acc, sequence_acc = line.strip().split()
                        if "|" in acc:
                            acc = acc.split("|")[1]
                        if acc not in sequences:
                            sequences[acc] = sequence_acc
                        else:
                            sequences[acc] += sequence_acc
                        if acc_human in acc:
                            acc_human = acc
                if not sequences.keys():
                    continue
                res = re.search(pattern, sequences[acc_human])
                begin = res.start()
                end = res.end()
                sequence = sequences[acc_human][begin:end]
                sequences = {i: j[begin:end] for i, j in sequences.items()}
                rev_seq = {}
                rev_seq_diffrence = {}
                print(file, sequences)

                for acc, seq in sequences.items():
                    if seq not in rev_seq:
                        rev_seq[seq] = set()
                    rev_seq[seq].add(acc)
                    rev_seq_diffrence[seq] = 0
                    for i in range(len(sequence)):
                        print("seq", acc, file, seq, i, sequence)
                        if seq[i] != sequence[i] and "canonical" not in rev_seq[seq]:
                            rev_seq_diffrence[seq] += 1
                can_seq = sequence
                for seq, acc in rev_seq.items():
                    if "canonical" in acc:
                        can_seq = seq
                position_H_by_protein = [5, 12, 13]
                first_pos = 0
                for i in range(len(can_seq)):
                    if first_pos < len(position_H_by_protein):
                        if i == position_H_by_protein[first_pos]:
                            first_pos += 1
                        elif i < position_H_by_protein[first_pos]:
                            if seq[i] == "-":
                                for first in range(first_pos, len(position_H_by_protein)):
                                    position_H_by_protein[first] += 1
                if len(rev_seq) == 2:
                    rev_seq = {i: j for i, j in rev_seq.items() if f"{j}" != "{'canonical_seq'}"}
                    saved = False
                    for seq, acc in rev_seq.items():
                        if "canonical" in f"{acc}":
                            saved = True
                            final_file.write(
                                f"{file.split('.')[0]};{seq};{acc};{seq[position_H_by_protein[0]] == 'H'};{seq[position_H_by_protein[1]] == 'H'};{seq[position_H_by_protein[2]] == 'H'}\n")
                    if not saved:
                        for seq, acc in rev_seq.items():
                            # if "canonical" not in f"{acc}":
                            if "canonical" not in f"{acc}":
                                saved = True
                                final_file.write(
                                    f"{file.split('.')[0]};{seq};{acc};{seq[position_H_by_protein[0]] == 'H'};{seq[position_H_by_protein[1]] == 'H'};{seq[position_H_by_protein[2]] == 'H'}\n")
                elif len(rev_seq) > 2:
                    max_diff = None
                    max_diff_seq = None
                    for seq, acc in rev_seq.items():
                        if len(acc) > 1 and "canonical" in f"{acc}":
                            max_diff_seq = seq
                            max_diff = 0
                            canonical_seq = seq
                        if "canonical" in f"{acc}":
                            canonical_seq = seq
                    for seq, acc in rev_seq.items():
                        if "canonical" not in f"{acc}":
                            if max_diff is None and max_diff_seq is None:
                                max_diff = distance(seq, canonical_seq)
                                max_diff_seq = seq
                            elif max_diff >= distance(seq, canonical_seq):
                                max_diff = distance(seq, canonical_seq)
                                max_diff_seq = seq
                    for seq, acc in rev_seq.items():
                        if max_diff_seq == seq:
                            final_file.write(
                                f"{file.split('.')[0]};{seq};{acc};{seq[position_H_by_protein[0]] == 'H'};{seq[position_H_by_protein[1]] == 'H'};{seq[position_H_by_protein[2]] == 'H'}\n")
                        o.write(f"{file.split('.')[0]};{seq};{acc};{rev_seq_diffrence[seq]};{max_diff_seq == seq}\n")
                elif len(rev_seq) == 1:
                    tmp_seq = list(rev_seq.keys())[0]
                    acc = rev_seq[sequence]
                    final_file.write(
                        f"{file.split('.')[0]};{tmp_seq};{acc};{seq[position_H_by_protein[0]] == 'H'};{seq[position_H_by_protein[1]] == 'H'};{seq[position_H_by_protein[2]] == 'H'}\n")
                with open(folder + "encoded_" + file, "w") as f:
                    for seq, acc in rev_seq.items():
                        f.write(f"{set_id}\t{seq}\n")
                        index_accs[set_id] = acc
                        set_id += 1

    with open(folder + "index.csv", "w") as f:
        for id_acc, accs in index_accs.items():
            f.write(f"{id_acc}\t{','.join(list(accs))}\n")


def align_final_abeta(file_in, file_out):
    with open(file_in) as f:
        with open("./tmp_file.fasta", "w") as f2:
            for line in f:
                if line.strip():
                    line_split = line.split(";")
                    f2.write(f">{line_split[0]}\n")
                    f2.write(f"{line_split[1]}\n")
    run_mafft("./tmp_file.fasta", "./tmp_file_mafft.aln")
    res = {}
    with open("./tmp_file_mafft.aln") as f:
        for line in f:
            if line.strip() and "CLUSTAL" not in line and len(line.split()) == 2:
                acc, ali = line.strip().split()
                if acc not in res:
                    res[acc] = ali
                else:
                    res[acc] += ali
    with open(file_out, "w") as f:
        for acc, seq in res.items():
            f.write(f"{acc};{seq}\n")


def get_ensemble(id):
    try:
        url = f"https://rest.ensembl.org/lookup/id/{id}?content-type=application/json;expand=1"
        res = requests.get(url)
        result = set([i["display_name"].split("-")[0] for i in res.json()["Transcript"]])
        if len(result) == 1:
            return list(result)[0].lower()
    except Exception:
        return


def get_ncbi(accessions):
    url = f" https://api.ncbi.nlm.nih.gov/datasets/v2/gene/accession/{accessions}"
    print(url)
    req = requests.get(url)
    resp_json = req.json()
    if "reports" in resp_json:
        for i in resp_json["reports"]:
            gene = i["gene"]
            print(gene["symbol"])
            print(f"https://www.ncbi.nlm.nih.gov/gene/{gene['gene_id']}")
            print(gene["description"])
            if "ensembl_gene_ids" in gene:
                ens = get_ensemble(gene["ensembl_gene_ids"][0])
                print(gene["ensembl_gene_ids"])
                if gene["symbol"] in ["app", "appa", "appb"]:
                    if ens == gene["symbol"]:
                        return gene["symbol"]
                else:
                    if ens in ["app", "appa", "appb"]:
                        return ens


def get_genes(source_file, save_file):
    with open(source_file) as f:
        with open(save_file, "w") as f2:
            for line in f:
                if line.startswith(">sp") or line.startswith(">tr"):
                    test = False
                    protein = line.split("|")[1]
                    url = f"https://rest.uniprot.org/uniprotkb/{protein}.xml"
                    print(url)
                    xml_file = requests.get(url)
                    url_text = xml_file.text.encode("utf-8")
                    soup = BeautifulSoup(url_text)
                    if line.split("|")[2] in ["app", "appa", "appb"]:
                        test = True
                    for gene in soup.find_all('gene'):
                        if gene.text.lower().strip() in ["appa", "appb"]:
                            if line.split("|")[2] not in ["app", "appa", "appb"]:
                                new_header = line.split("|")[:2] + [gene.text.lower().strip()] + line.split("|")[2:]
                                new_header = "|".join(new_header)
                                line = new_header
                                test = True

                    if not test:
                        ensembl = soup.find_all('dbreference', type="Ensembl")
                        for ens in ensembl:
                            properties = ens.find_all("property", type="gene ID")
                            for prop in properties:
                                ensembl_id = prop["value"]
                                gene_name = get_ensemble(ensembl_id.split(".")[0])
                                if gene_name in ["app", "appa", "appb"]:
                                    new_header = line.split("|")[:2] + [gene_name] + line.split("|")[2:]
                                    new_header = "|".join(new_header)
                                    line = new_header

                        for embl in soup.find_all('dbreference', type="EMBL"):
                            embl_id = embl["id"]
                            print(f"https://www.ebi.ac.uk/ena/browser/view/{embl_id}")
                            gene_name = input().strip()
                            if gene_name:
                                new_header = line.split("|")[:2] + [gene_name] + line.split("|")[2:]
                                new_header = "|".join(new_header)
                                line = new_header
                                print(new_header)
                        for refseq in soup.find_all('dbreference', type="RefSeq"):
                            print(refseq)
                            gene_name = input().strip()
                            if gene_name:
                                new_header = line.split("|")[:2] + [gene_name] + line.split("|")[2:]
                                new_header = "|".join(new_header)
                                line = new_header
                elif line.startswith(">"):
                    refseq_id = line[1:].split()[0]
                    if "_" in line:
                        print(line.split("_")[-1].split()[0])
                        if line.split("_")[-1].split()[0] in ["app", "appa", "appb"]:
                            continue
                        else:
                            print(f"https://www.ncbi.nlm.nih.gov/protein/{refseq_id}")
                            gene_name = get_ncbi(refseq_id)
                            if gene_name:
                                new_header = [line.split()[0] + f"_{gene_name}"] + line.split()[1:]
                                new_header = " ".join(new_header)
                                line = new_header
                    else:
                        print(f"https://www.ncbi.nlm.nih.gov/protein/{refseq_id}")
                        gene_name = get_ncbi(refseq_id)
                        if gene_name:
                            new_header = [line.split()[0] + f"_{gene_name}"] + line.split()[1:]
                            new_header = " ".join(new_header)
                            line = new_header
                f2.write(line)


def get_actino_sequences(fasta_sequences_by_headers, path):
    ncbi = NCBITaxa()
    actino_fasta = {}
    for acc, data_fasta in fasta_sequences_by_headers.items():
        if ">sp" in data_fasta["header"] or ">tr" in data_fasta["header"]:
            tax = get_tax_uniprot(acc, data_fasta["header"])
        else:
            organism = data_fasta["header"].split("[")[-1].split("]")[0]
            tax = get_organism_by_refseq(organism, data_fasta)
        print(acc, tax)
        if tax is None:
            print(data_fasta)
            print("tax_id:")
            tax = input()
        if 7898 in ncbi.get_lineage(tax):
            actino_fasta[acc] = data_fasta
    with open(path, "w") as f:
        for acc, seq_data in actino_fasta.items():
            f.write(seq_data["header"])
            f.write(seq_data["seq"])


def read_files_with_sequences(file):
    res = {}
    new_str = ""
    name = ""
    header = ""
    with open(file) as f:
        for line in f.readlines():
            if line.startswith(">"):
                if header:
                    res[name] = {"header": header, "seq": new_str}
                    if ">sp" in line or ">tr" in line:
                        name = line.split("|")[1]
                    else:
                        name = line.split()[0][1:].replace("_appa", "").replace("_appb", "").replace("_app", "")
                header = line
                new_str = ""
            else:
                new_str += line
    if name:
        res[name] = {"header": header, "seq": new_str}
    return res


def tuning_cd_hit(file):
    c_list = [0.7, 0.8, 0.9, 0.98]
    s_list = [0.7, 0.8, 0.9]
    aL_list = [0.7, 0.8, 0.9]
    for c in c_list:
        for s in s_list:
            for aL in aL_list:
                command = f"cd-hit -c {c} -n 5 -G 1 -i {file}  -o  tmp_c_{c}_s_{s}_aL_{aL} -s {s} -aL {aL} -p 1 -T 0 -g 1"
                os.system(command)


def change_gene(old_gene, new_gene, cl_data, fasta):
    for acc in cl_data[old_gene]:
        if ">sp" in fasta[acc]["header"] or ">tr" in fasta[acc]["header"]:
            if f"|{new_gene}|" not in fasta[acc]["header"]:
                fasta[acc]["header"] = fasta[acc]["header"].replace(f"{acc}|", f"{acc}|{new_gene}|")
        else:
            if f"_{new_gene}" not in fasta[acc]["header"]:
                fasta[acc]["header"] = fasta[acc]["header"].replace(f"{acc}", f"{acc}_{new_gene}")
    for acc in cl_data["app"]:
        if ">sp" in fasta[acc]["header"] or ">tr" in fasta[acc]["header"]:
            if f"|{new_gene}|" not in fasta[acc]["header"]:
                fasta[acc]["header"] = fasta[acc]["header"].replace(f"{acc}|", f"{acc}|{new_gene}|")
        else:
            if f"_{new_gene}" not in fasta[acc]["header"]:
                fasta[acc]["header"] = fasta[acc]["header"].replace(f"{acc}", f"{acc}_{new_gene}")
    for acc in cl_data["lack"]:
        if ">sp" in fasta[acc]["header"] or ">tr" in fasta[acc]["header"]:
            if f"|{new_gene}|" not in fasta[acc]["header"]:
                fasta[acc]["header"] = fasta[acc]["header"].replace(f"{acc}|", f"{acc}|{new_gene}|")
        else:
            if f"_{new_gene}" not in fasta[acc]["header"]:
                fasta[acc]["header"] = fasta[acc]["header"].replace(f"{acc}", f"{acc}_{new_gene}")
    return fasta


def change_genes_names(last_file, fasta, file_proteins_appa_appb_with_genes):
    c_list = [0.98, 0.95, 0.9, 0.8, 0.85]
    s_list = [0.98, 0.95, 0.9, 0.8, 0.85]
    aL_list = [0.98, 0.95, 0.9, 0.8, 0.85]
    for c in c_list:
        for s in s_list:
            for aL in aL_list:
                command = f"cd-hit -c {c} -n 5 -G 1 -i {last_file}  -o  tmp_c_{c}_s_{s}_aL_{aL} -s {s} -aL {aL} -p 1 -T 0 -g 1"
                os.system(command)
                result = {}
                cluster_name = None
                new_data = {"lack": set(), "app": set(), "appa": set(), "appb": set()}
                # with open(f"./new_tmp/new_tmp_c_{c}_s_{s}_aL_{aL}.clstr") as f:
                with open(f"./tmp_c_{c}_s_{s}_aL_{aL}.clstr") as f:
                    for line in f:
                        if line.startswith(">"):
                            if cluster_name:
                                result[cluster_name] = new_data
                            cluster_name = line.strip()
                            new_data = {"lack": set(), "app": set(), "appa": set(), "appb": set()}
                        else:
                            gene = "lack"
                            line_splited = line.split(" ")
                            acc_data = line_splited[1]
                            # print(f"tmp_c_{c}_s_{s}_aL_{aL}.clstr", cluster_name, acc_data)

                            if "tr" in acc_data or "sp" in acc_data:
                                acc = acc_data.split("|")[1]
                                gene = acc_data.split("|")[2]
                                # print(acc_data, gene)

                            else:
                                gene = "lack"
                                if "app" in acc_data or "appa" in acc_data or "appb" in acc_data:
                                    acc = acc_data[1:].rsplit("_", 1)[0].replace("_appa", "").replace("_appb",
                                                                                                      "").replace(
                                        "_app", "").replace("?", "")
                                    gene = acc_data[1:].rsplit("_")[-1][:-3]
                                    # print(acc_data, gene)

                                else:
                                    acc = acc_data[1:-3]
                            if gene not in ["app", "appa", "appb"]:
                                gene = "lack"
                            new_data[gene].add(acc)
                # print(f"tmp_c_{c}_s_{s}_aL_{aL}.clstr", result)
                for cl_name, cl_data in result.items():
                    if cl_data["appa"] and cl_data["appb"]:
                        if len(cl_data["appb"]) < len(cl_data["appa"]):
                            print(f"tmp_c_{c}_s_{s}_aL_{aL}.clstr", "mistake!", cl_name, "b:", cl_data["appb"], "a:",
                                  cl_data["appa"])
                        else:
                            print(f"tmp_c_{c}_s_{s}_aL_{aL}.clstr", "mistake!", cl_name, "a:", cl_data["appa"], "b:",
                                  cl_data["appb"])
                    if len(cl_data["appa"]) + len(cl_data["appb"]) + len(cl_data["lack"]) + len(cl_data["app"]) > 0:
                        if len(cl_data["appa"]) / (
                                len(cl_data["appa"]) + len(cl_data["appb"]) + len(cl_data["lack"]) + len(
                            cl_data["app"])) > 0.5 and len(cl_data["appa"]) > len(cl_data["appb"]):
                            fasta = change_gene(old_gene="appb", new_gene="appa", cl_data=cl_data, fasta=fasta)
                        elif len(cl_data["appb"]) / (
                                len(cl_data["appa"]) + len(cl_data["appb"]) + len(cl_data["lack"]) + len(
                            cl_data["app"])) > 0.5 and len(cl_data["appb"]) > len(cl_data["appa"]):
                            fasta = change_gene(old_gene="appa", new_gene="appb", cl_data=cl_data, fasta=fasta)
                all = 0
                with_gene = 0
                with open(file_proteins_appa_appb_with_genes, "w") as f2:
                    for acc, fasta_data in fasta.items():
                        all += 1
                        if "appa" in fasta_data["header"] or "appb" in fasta_data["header"]:
                            with_gene += 1
                        else:
                            f2.write(fasta_data["header"])
                            f2.write(fasta_data["seq"])
    return fasta


def divide_appa_appb_by_files(fasta, appa_folder, appb_folder, canonical_sequence):
    info_tax = {}
    for acc, data_fasta in fasta.items():
        if ">sp" in data_fasta["header"] or ">tr" in data_fasta["header"]:
            tax = get_tax_uniprot(acc)
            gene = data_fasta["header"].split("|")[2]
        else:
            organism = data_fasta["header"].split("[")[-1].split("]")[0]
            tax = get_organism_by_refseq(organism)
            gene = "ap" + data_fasta["header"].split(" ")[0].split("_ap", 1)[-1]
            if "appa" in gene:
                gene = "appa"
            if "appb" in gene:
                gene = "appb"
        if tax not in info_tax:
            info_tax[tax] = {"appa": set(), "appb": set(), "app": set(), "other": set()}
        print(gene)
        if gene in ["appa", "appb", "app"]:
            info_tax[tax][gene].add(acc)
        else:
            info_tax[tax]["other"].add(acc)
        gene = None
    with open("summary_actino.csv", "w") as f:
        for tax, info in info_tax.items():
            f.write(
                f"{tax};{','.join(list(info['appa']))};{','.join(list(info['appb']))};{','.join(list(info['app']))};{','.join(list(info['other']))}\n")
    # tax;appa;appb;app;other
    for tax, info in info_tax.items():
        for prot in info['appa']:
            if not os.path.exists(os.path.join(appa_folder, f"{tax}.fasta")):
                with open(os.path.join(appa_folder, f"{tax}.fasta"), "w") as f2:
                    f2.write(canonical_sequence)
            with open(os.path.join(appa_folder, f"{tax}.fasta"), "a") as f:
                f.write(fasta[prot]["header"])
                f.write(fasta[prot]["seq"])
        for prot in info['appb']:
            if not os.path.exists(os.path.join(appb_folder, f"{tax}.fasta")):
                with open(os.path.join(appb_folder, f"{tax}.fasta"), "w") as f2:
                    f2.write(canonical_sequence)
            with open(os.path.join(appb_folder, f"{tax}.fasta"), "a") as f:
                f.write(fasta[prot]["header"])
                f.write(fasta[prot]["seq"])


@click.command()
@click.option('--fasta_uniprot_file', default="./data/abeta_blast_uniprot.fasta", help='Path to UniProt fasta file.')
@click.option('--fasta_refseq_file', default="./data/abeta_blast_refseq.fasta", help='Path to RefSeq fasta file.')
@click.option('--result_folder', default="./result/", help='Path to result folder.')
@click.option('--canonical_sequence', default="./data/canonical_sequence.fasta",
              help='Path to file with reference_sequence.')
def run_blast(fasta_uniprot_file, fasta_refseq_file, result_folder, canonical_sequence):
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    canonical_sequence = "\n".join([i for i in open(canonical_sequence).read().split("\n") if i.strip()])
    file_refseq_alignment = os.path.join(result_folder, "alignment_refseq.aln")
    file_refseq_proper_alignment = os.path.join(result_folder, "alignment_refseq_AB.aln")
    file_refseq_wrong_alignment = os.path.join(result_folder, "alignment_refseq_AB_excluded.aln")
    file_refseq_organisms = os.path.join(result_folder, "organisms_refseq.csv")

    files_uniprot_alignment = os.path.join(result_folder, "alignment_uniprot.aln")
    file_uniprot_proper_alignment = os.path.join(result_folder, "alignment_uniprot_AB.aln")
    file_uniprot_wrong_alignment = os.path.join(result_folder, "alignment_uniprot_AB_excluded.aln")
    folder_organisms = os.path.join(result_folder, "organism")
    if not os.path.exists(folder_organisms):
        os.mkdir(folder_organisms)
    folder_organisms_with_reference = os.path.join(result_folder, "organism_updated")
    if not os.path.exists(folder_organisms_with_reference):
        os.mkdir(folder_organisms_with_reference)

    organims_with_isoforms = os.path.join(result_folder, "problematic_organism2.csv")
    result_file = os.path.join(result_folder, "final_file.csv")
    result_file_aligned = os.path.join(result_folder, "final_file_aligned.csv")
    ox_sets, fasta_sequences = dict(), dict()

    appa_appb_fasta_file = os.path.join(result_folder, "appa_appb.fasta")
    appa_appb_gene_fasta_file = os.path.join(result_folder, "appa_appb_genes.fasta")
    appa_appb_file_proteins_with_genes = os.path.join(result_folder, '"appa_appb_only_genes.fasta"')

    appa_organisms_folder = os.path.join(result_folder, "appa_organisms")
    appb_organisms_folder = os.path.join(result_folder, "appb_organisms")

    if not os.path.exists(appa_organisms_folder):
        os.mkdir(appa_organisms_folder)
    if not os.path.exists(appb_organisms_folder):
        os.mkdir(appb_organisms_folder)

    if fasta_refseq_file is not None:
        run_mafft(fasta_refseq_file, file_refseq_alignment)
        search_APP_localisation(file_aln=file_refseq_alignment,
                                file_out_aln=file_refseq_proper_alignment,
                                file_out_aln_excluded=file_refseq_wrong_alignment)
        ox_sets, fasta_sequences = get_organisms_refseq(file_aln_refseq=file_refseq_proper_alignment,
                                                        fasta_file_refseq=fasta_refseq_file,
                                                        organism_file=file_refseq_organisms, ox_sets=ox_sets,
                                                        fasta_sequences=fasta_sequences)
    if fasta_uniprot_file is not None:
        run_mafft(fasta_uniprot_file, files_uniprot_alignment)
        search_APP_localisation(file_aln=files_uniprot_alignment,
                                file_out_aln=file_uniprot_proper_alignment,
                                file_out_aln_excluded=file_uniprot_wrong_alignment)
        ox_sets, fasta_sequences = get_organisms_uniprot(file_aln_uniprot=file_uniprot_proper_alignment,
                                                         fasta_file_uniprot=fasta_uniprot_file,
                                                         ox_sets=ox_sets, fasta_sequences=fasta_sequences)
    fasta_sequences_by_headers = {i: {"header": j.split("\n")[0], "seq": j.split("\n")[0]} for i, j in fasta_sequences}
    get_actino_sequences(fasta_sequences_by_headers, appa_appb_fasta_file)

    get_genes(appa_appb_fasta_file, appa_appb_gene_fasta_file)
    appa_appb_sequences = read_files_with_sequences(appa_appb_gene_fasta_file)
    tuning_cd_hit(appa_appb_gene_fasta_file)
    changed_genes_fasta = change_genes_names(appa_appb_gene_fasta_file, appa_appb_sequences,
                                             appa_appb_file_proteins_with_genes)

    divide_appa_appb_by_files(changed_genes_fasta, appa_organisms_folder, appb_organisms_folder, canonical_sequence)

    make_mafft_per_organism(appb_organisms_folder)
    make_mafft_per_organism(appa_organisms_folder)
    encode_mafft_find_amyloid_per_organism(appa_organisms_folder, os.path.join(appa_organisms_folder, result_file),
                                           organims_with_isoforms)
    encode_mafft_find_amyloid_per_organism(appb_organisms_folder, os.path.join(appb_organisms_folder, result_file),
                                           organims_with_isoforms)

    divide_by_organisms(folder_organisms, ox_sets, fasta_sequences)
    update_organisms(fasta_sequences=fasta_sequences, ox_sets=ox_sets, out_folder=folder_organisms_with_reference,
                     canonical_sequence=canonical_sequence)
    make_mafft_per_organism(folder_organisms_with_reference)
    encode_mafft_find_amyloid_per_organism(folder_organisms_with_reference, result_file, organims_with_isoforms)
    file_names = os.listdir(appb_organisms_folder)
    for file_name in file_names:
        shutil.move(os.path.join(appb_organisms_folder, file_name), folder_organisms_with_reference)
    align_final_abeta(result_file, result_file_aligned)


if __name__ == "__main__":
    run_blast()
