import requests


def get_sequence(unirefacc):
    url = f"https://rest.uniprot.org/uniref/{unirefacc}.fasta"
    req = requests.get(url)
    if req.status_code == 200:
        return req.text
    print(url)


def read_blast(file):
    with open(file) as f:
        for line in f:
            if line.startswith(">"):
                unirefacc = line.split()[0][1:]
                yield get_sequence(unirefacc)


def run():
    file = "C:\\Users\\IBB\\PycharmProjects\\abeta_evolution\\data2\\nr_db_blast\\abeta_blast_uniprot_locally.txt"
    file_out = "./abeta_blast_uniprot.fasta"
    blast_res = read_blast(file)

    with open(file_out, "w") as f:
        for protein in blast_res:
            # print(protein)
            if protein:
                f.write(protein.replace("UniRef90_", ""))


if __name__ == "__main__":
    run()
