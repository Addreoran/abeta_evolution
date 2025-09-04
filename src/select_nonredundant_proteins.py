import os
import re
import requests
from Levenshtein import distance
import click


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
        for l in f:
            if l.strip():
                acc, tax_id = l.strip().split(";")
                result[acc] = tax_id
    return result


def get_organism_by_refseq(refseq_id):
    from ete3 import NCBITaxa
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
    for l in req.text.split("\n"):
        if l.startswith(">"):
            return l.split("OX=")[-1].split()[0]


def get_organisms_uniprot(file_aln_uniprot, fasta_file_uniprot, fasta_sequences=None, ox_sets=None):
    acc_ox = {}
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
                    except:
                        print(f"TaxId not found. Give me TaxId of {acc}:")
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


def update_organisms(fasta_sequences, ox_sets, out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    for ox, accs in ox_sets.items():
        with open(f"{out_folder}/{ox}.fasta", "w") as f:
            f.write(">1|canonical_seq\nMLPGLALLLLAAWTARALEVPTDGNAGLLAEPQIAMFCGRLNMHMNVQNGKWDSDPSGTK\n"
                    "TCIDTKEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPHFVIPYRCLVG\n"
                    "EFVSDALLVPDKCKFLHQERMDVCETHLHWHTVAKETCSEKSTNLHDYGMLLPCGIDKFR\n"
                    "GVEFVCCPLAEESDNVDSADAEEDDSDVWWGGADTDYADGSEDKVVEVAEEEEVAEVEEE\n"
                    "EADDDEDDEDGDEVEEEAEEPYEEATERTTSIATTTTTTTESVEEVVREVCSEQAETGPC\n"
                    "RAMISRWYFDVTEGKCAPFFYGGCGGNRNNFDTEEYCMAVCGSAMSQSLLKTTQEPLARD\n"
                    "PVKLPTTAASTPDAVDKYLETPGDENEHAHFQKAKERLEAKHRERMSQVMREWEEAERQA\n"
                    "KNLPKADKKAVIQHFQEKVESLEQEAANERQQLVETHMARVEAMLNDRRRLALENYITAL\n"
                    "QAVPPRPRHVFNMLKKYVRAEQKDRQHTLKHFEHVRMVDPKKAAQIRSQVMTHLRVIYER\n"
                    "MNQSLSLLYNVPAVAEEIQDEVDELLQKEQNYSDDVLANMISEPRISYGNDALMPSLTET\n"
                    "KTTVELLPVNGEFSLDDLQPWHSFGADSVPANTENEVEPVDARPAADRGLTTRPGSGLTN\n"
                    "IKTEEISEVKMDAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIATVIVITL\n"
                    "VMLKKKQYTSIHHGVVEVDAAVTPEERHLSKMQQNGYENPTYKFFEQMQN\n")
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
        no_duble = 0
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
            for l in f:
                if l.strip():
                    line = l.split(";")
                    f2.write(f">{line[0]}\n")
                    f2.write(f"{line[1]}\n")
    run_mafft("./tmp_file.fasta", "./tmp_file_mafft.aln")
    res = {}
    with open("./tmp_file_mafft.aln") as f:
        for l in f:
            if l.strip() and "CLUSTAL" not in l and len(l.split()) == 2:
                print(l, l.strip().split())
                acc, ali = l.strip().split()
                if acc not in res:
                    res[acc] = ali
                else:
                    res[acc] += ali
    with open(file_out, "w") as f:
        for acc, seq in res.items():
            f.write(f"{acc};{seq}\n")


@click.command()
@click.option('--fasta_uniprot_file', default="./data/abeta_blast_uniprot.fasta", help='Path to UniProt fasta file.')
@click.option('--fasta_refseq_file', default="./data/abeta_blast_refseq.fasta", help='Path to RefSeq fasta file.')
@click.option('--result_folder', default="./result/", help='Path to result folder.')
def run_blast(fasta_uniprot_file, fasta_refseq_file, result_folder):
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    file_refseq_alignment = os.path.join(result_folder, "alignment_refseq.aln")
    file_refseq_proper_alignment = os.path.join(result_folder, "alignment_refseq_AB.aln")
    file_refseq_wrong_alignment = os.path.join(result_folder, "alignment_refseq_AB_excluded.aln")
    file_refseq_organisms = os.path.join(result_folder, "organisms_refseq.csv")

    files_uniprot_alignment = os.path.join(result_folder, "alignment_uniprot.aln")
    file_uniprot_proper_alignment = os.path.join(result_folder, "alignment_uniprot_AB.aln")
    file_uniprot_wrong_alignment = os.path.join(result_folder, "alignment_uniprot_AB_excluded.aln")

    folder_organisms = os.path.join(result_folder, "organism")
    folder_organisms_with_reference = os.path.join(result_folder, "organism_updated/")
    organims_with_isoforms = os.path.join(result_folder, "problematic_organism2.csv")
    result_file = os.path.join(result_folder, "final_file.csv")
    result_file_aligned = os.path.join(result_folder, "final_file_aligned.csv")
    ox_sets, fasta_sequences = dict(), dict()

    if fasta_refseq_file is not None:
        run_mafft(fasta_refseq_file, file_refseq_alignment)
        search_APP_localisation(file_aln=file_refseq_alignment,
                                file_out_aln=file_refseq_proper_alignment,
                                file_out_aln_excluded=file_refseq_wrong_alignment)
        ox_sets, fasta_sequences = get_organisms_refseq(file_aln_refseq=file_refseq_proper_alignment,
                                                        fasta_file_refseq=fasta_refseq_file,
                                                        organism_file=file_refseq_organisms, ox_sets=ox_sets,
                                                        fasta_sequences=fasta_sequences)
    # todo: add chondroichties analyse
    if fasta_uniprot_file is not None:
        run_mafft(fasta_uniprot_file, files_uniprot_alignment)
        search_APP_localisation(file_aln=files_uniprot_alignment,
                                file_out_aln=file_uniprot_proper_alignment,
                                file_out_aln_excluded=file_uniprot_wrong_alignment)
        ox_sets, fasta_sequences = get_organisms_uniprot(file_aln_uniprot=file_uniprot_proper_alignment,
                                                         fasta_file_uniprot=fasta_uniprot_file,
                                                         ox_sets=ox_sets, fasta_sequences=fasta_sequences)
    # todo: add chondroichties analyse
    divide_by_organisms(folder_organisms, ox_sets, fasta_sequences)
    update_organisms(fasta_sequences=fasta_sequences, ox_sets=ox_sets, out_folder=folder_organisms_with_reference)
    make_mafft_per_organism(folder_organisms_with_reference)
    encode_mafft_find_amyloid_per_organism(folder_organisms_with_reference, result_file, organims_with_isoforms)
    align_final_abeta(result_file, result_file_aligned)


if __name__ == "__main__":
    run_blast()
