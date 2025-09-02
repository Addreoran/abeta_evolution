import os.path

from ete3 import NCBITaxa

from protein import Protein
import click
import weblogo


def get_groups():
    groups = dict(
        fish=(7777, 7898, 118072),
        bony_fish=(7878, 118072,),
        mammals=(40674),
        aves=(8782),
        reptiles=(8504, 1294634, 8459),
        amphibians=(8292)
    )
    return groups


def read_file(file):
    file = open(file, "r")
    proteins = []
    for i in file.readlines():
        print(i)
        protein = Protein(uniprot_id=i.strip().split(";")[-1].replace(".aln", ""),
                          organism_id=i.strip().split(";")[0].replace('\n', '').replace(".aln", ""),
                          seq=i.split(";")[1])
        protein.line = i
        proteins.append(protein)
    return proteins


def get_proteins_from_group(proteins, group):
    ncbi = NCBITaxa()
    proteins2 = [i for i in proteins if i.get_organism(ncbi, group) in group]
    return proteins2


def get_group(proteins, taxId):
    ncbi = NCBITaxa()
    proteins_new = [i for i in proteins if i.get_organism(ncbi, typ="general", selected=taxId) in taxId]
    return proteins_new


def save_group(proteins, file):
    with open(file, "w") as f:
        for p in proteins:
            f.write(p.line)


def read_divided_file(file):
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
            if prot_id != "taxId":
                f.write(f">{prot_id}\n")
                f.write(f"{sequence.replace('-', '')[:28]}\n")


def run_mafft(file, fasta_file):
    command = f'"/usr/bin/mafft"  --anysymbol --auto --clustalout --reorder "{fasta_file}" > "{file}"'
    os.system(command)


def parse_mafft(file, file_out):
    sequences = {}
    with open(file) as f:
        # 1 otworzyć i spisać sekwencje
        for line in f.readlines():
            if line.strip() and "CLUSTAL" not in line:
                if len(line.strip().split()) == 2:
                    acc, sequence_acc = line.strip().split()
                    if acc not in sequences:
                        sequences[acc] = sequence_acc
                    else:
                        sequences[acc] += sequence_acc
    with open(file_out, "w") as f:
        for acc, seq in sequences.items():
            f.write(f"{acc};")
            f.write(f"{seq}\n")


def draw_logo(save_file, file):
    with open("logo.cln", "w") as f1:
        f1.write("CLUSTAL 2.1 multiple sequence alignment\n\n")
        sab = {}
        line_no = 0
        with open(file) as f:
            for l in f.readlines():
                if not l.startswith("taxId"):
                    if l.strip() and "*" not in l:
                        line = l.strip().split(";")
                        uniprot = line_no
                        line_no += 1
                        seq = line[1]
                        sab[uniprot] = seq[:16]
            positions = [0 for i in range(len(seq))]
            popular_positions = []
            for uniprot, seq in sab.items():
                for i in range(len(seq)):
                    if seq[i] != "-":
                        positions[i] += 1
                    else:
                        print(uniprot, i)
            for e, i in enumerate(positions):
                if i / len(sab.keys()) > 0.1:
                    popular_positions.append(e)
            for uniprot, seq in sab.items():
                seq_save = "".join([seq[i] for i in popular_positions])
                f1.write(f"{uniprot}\t{seq_save}\n")
    fin = open('logo.cln')
    seqs = weblogo.read_seq_data(fin)
    logodata = weblogo.LogoData.from_seqs(seqs)
    logooptions = weblogo.LogoOptions()
    logooptions.unit_name = "probability"
    logooptions.scale_width = False
    logoformat = weblogo.LogoFormat(logodata, logooptions)
    # logoformat.logo_width = logoformat.logo_width * 5
    # logoformat.logo_height = logoformat.logo_height * 5
    jpeg = weblogo.logo_formatter.svg_formatter(logodata, logoformat)
    with open(save_file, "ab") as a:
        a.write(jpeg)


def run_logo_steps(proteins, taxes, result_folder, name):
    proteins2 = get_group(proteins, taxId=taxes)
    save_group(proteins2, os.path.join(result_folder, f"peptides_{name}.csv"))

    res = read_divided_file(os.path.join(result_folder, f"peptides_{name}.csv"))
    save_as_fasta(os.path.join(result_folder, f"peptides_{name}.fasta"), res)
    run_mafft(os.path.join(result_folder, f"peptides_{name}.aln"),
              os.path.join(result_folder, f"peptides_{name}.csv"))
    parse_mafft(os.path.join(result_folder, f"peptides_{name}.aln"),
                os.path.join(result_folder, f"peptides_{name}_aligned.csv"))
    draw_logo(os.path.join(result_folder, f"peptides_{name}.svg"),
              os.path.join(result_folder, f"peptides_{name}_aligned.csv"))


@click.option('--result_file', default="./result/final_file_aligned.csv",
              help='Path to finale file from script "select_nonredundant_proteins.py".')
@click.option('--TaxIds', default=[],
              help='List of Taxonomy IDs to plot. If empty, it divide and plot default groups of animals.')
@click.option('--result_folder', default="./WebLogo/",
              help='Folder with results.')
def divide_organisms_by_groups(result_file, TaxIds, result_folder):
    proteins = read_file(result_file)
    if not TaxIds:
        groups = get_groups()
        for name, taxes in groups.items():
            run_logo_steps(proteins, taxes, result_folder, name)
    else:
        run_logo_steps(proteins, TaxIds, result_folder, '_'.join(TaxIds))
