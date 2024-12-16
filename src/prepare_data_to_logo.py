import os


def read_file(file):
    result = {}
    with open(file) as f:
        for l in f:
            if l.strip():
                line = l.strip().split(";")
                if line:
                    result[line[0]] = line[1]
    return result


def save_as_fasta(file, result):
    with open(file, "w") as f:
        for prot_id, sequence in result.items():
            f.write(f">{prot_id}\n")
            f.write(f"{sequence.replace('-', '')}\n")


#
def run_mafft(file, fasta_file):
    command = f'"/usr/bin/mafft"  --anysymbol --auto --clustalout --reorder "{fasta_file}" > "{file}"'
    os.system(command)


#
#
def parse_mafft(file, file_out):
    sequences = {}
    with open(file) as f:
        # 1 otworzyć i spisać sekwencje
        for line in f.readlines():
            if line.strip() and "CLUSTAL" not in line:
                if len(line.strip().split())==2:
                    acc, sequence_acc = line.strip().split()
                    if acc not in sequences:
                        sequences[acc] = sequence_acc
                    else:
                        sequences[acc] += sequence_acc
    with open(file_out, "w") as f:
        for acc, seq in sequences.items():
            f.write(f"{acc};")
            f.write(f"{seq}\n")


files = ["../data/final_file.csv",
         "../data/final_results/divided_groups_without_gene_redundancy_fish_Q.csv",
         f"../data/final_results/divided_groups_without_gene_redundancy_Chondrichtlhyes.csv",
         "../data/final_results/divided_groups_without_gene_redundancy_Birds.csv",
         f"../data/final_results/divided_groups_without_gene_redundancy_Amphibians.csv",
         f"../data/final_results/divided_groups_without_gene_redundancy_reptiles.csv",
         "../data/final_results/divided_groups_without_gene_redundancy_fish.csv",
         "../data/final_results/divided_groups_without_gene_redundancy_mammals.csv",
         "../data/final_results/divided_groups_without_gene_redundancy_Actinopterygii.csv"]
files_out = ["../data/final_file.fasta",
             "../data/final_results/divided_groups_without_gene_redundancy_fish_Q.fasta",
             f"../data/final_results/divided_groups_without_gene_redundancy_Chondrichtlhyes.fasta",
             "../data/final_results/divided_groups_without_gene_redundancy_Birds.fasta",
             f"../data/final_results/divided_groups_without_gene_redundancy_Amphibians.fasta",
             f"../data/final_results/divided_groups_without_gene_redundancy_reptiles.fasta",
             "../data/final_results/divided_groups_without_gene_redundancy_fish.fasta",
             "../data/final_results/divided_groups_without_gene_redundancy_mammals.fasta",
             "../data/final_results/divided_groups_without_gene_redundancy_Actinopterygii.fasta"]
for i in range(len(files)):
    res = read_file(files[i])
    save_as_fasta(files_out[i], res)
    run_mafft(files_out[i].rsplit(".", 1)[0] + ".aln", files_out[i])
    parse_mafft(files_out[i].rsplit(".", 1)[0] + ".aln", files_out[i].rsplit(".", 1)[0] + "_aligned.csv")
