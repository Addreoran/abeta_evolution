import os
import re

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
    os.system(f'"/usr/bin/mafft"  --auto --clustalout --reorder "{file}" > "{out}"')


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
    sequences_included = {i: j for i, j in sequences.items() if j[begin:end].count("-") > 5}
    sequences_excluded = {i: j for i, j in sequences.items() if j[begin:end].count("-") > 5}

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


if __name__ == "__main__":
    # 1) Pobranie białek powstających z genu APP
    # https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28gene%3AAPP%29+AND+%28protein_name%3AAmyloid-beta%29%29
    # 2) Usunięcie białek, które nie powstają z APP
    get_APP_protein(file_uniprot="../data/uniprot_APP.fasta", file_app="../data/uniprot_only_APP.fasta")
    # 3) szukamy fragmentów abeta
    # 3A) mafft
    run_mafft("../data/uniprot_only_APP.fasta", "../data/alignment_APP.aln")
    # 3b) lokalizacja abety i zapisanie do pliku
    sequences = search_APP_localisation(file_aln="../data/alignment_APP.aln", file_out_aln="../data/alignment_AB.aln",
                                        file_out_aln_excluded="../data/alignment_AB_excluded.aln")
    # 3c) zapisanie fasta do pliku
    aln_to_fasta(file_fasta_all="../data/uniprot_only_APP.fasta",
                 file_fasta="../data/uniprot_AB.fasta",
                 sequences=sequences)
    # 4) jackhmmer

    # 4a) pobranie bazy trembl+sp

    # 4b) puszczenie jackhmmer z input jako "../data/uniprot_AB.fasta"

    # 5) usunięcie redundancji
