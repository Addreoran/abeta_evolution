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
                    organism = line.split("OX=")[-1].split("OS=")[-1].split()[0]
                    acc_headers[acc] = line
    seq_data = {}
    with open(file_aln) as f:
        for line in f.readlines():
            if line.strip():
                line = line.split()
                acc = line[0].split("|")[1]
                sequence = line[1].replace("-", "")
                seq_data[acc] = sequence
    with open(fasta_ab, "w") as f:
        for acc, seq in seq_data.items():
            f.write(f">{acc}\n")
            f.write(f"{seq}\n")

    pass


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
    # 3d) sprawdzić czy wykluczenia są ok
    # - pobrać nazwy białek i orgainzmy,
    # jeśli takie samo już jest w moim zestawie białek,
    # to nie wyrzucam, jeśli nie ma, to sprawdzam czy na pewno powinno zostać wyrzucone
    get_excluded_names(file_fasta="../data/uniprot_only_APP.fasta",
                       file_excluded="../data/alignment_AB_excluded.aln",
                       file_aln="../data/alignment_AB.aln",
                       analise_file="../data/excluded_acc_analyse.csv")
    # 4) jackhmmer
    # 4a) pobranie bazy trembl+sp
    # get_uniprot("../data/uniprot.fasta")
    # 4b) pobranie fasta z align
    get_fasta_of_aln(file_aln="../data/alignment_AB.aln",
              fasta_all="../data/uniprot_APP.fasta",
              fasta_ab="../data/AB.fasta")
    # 4b) puszczenie jackhmmer z input jako "../data/uniprot_AB.fasta"

    # 5) mafft po raz drugi z usunięciem braków w alignmencie lub inne białka

    # 6) usunięcie redundancji
