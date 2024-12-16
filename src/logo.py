# from weblogo import logo_formatter, read_seq_data, LogoData, LogoOptions, LogoFormat
import weblogo


def draw_logo(save_file, file):
    with open("logo.cln", "w") as f1:
        f1.write("CLUSTAL 2.1 multiple sequence alignment\n\n")
        sab = {}
        line_no = 0
        with open(file) as f:
            for l in f.readlines():
                line = l.strip().split(";")
                uniprot = line_no
                line_no += 1
                seq = line[1]
                sab[uniprot] = seq
                # print(len(seq))
                # print(seq)
            positions = [0 for i in range(len(seq))]
            popular_positions = []
            for uniprot, seq in sab.items():
                for i in range(len(seq)):
                    print(uniprot, i, len(seq), len(positions), seq[i], seq, positions)
                    if seq[i] != "-":
                        positions[i] += 1
                    else:
                        print(uniprot, i)
            for e, i in enumerate(positions):
                if i / len(sab.keys()) > 0.5:
                    popular_positions.append(e)
                print(i, len(sab.keys()), i / len(sab.keys()), popular_positions)
            print(popular_positions)
            for uniprot, seq in sab.items():
                seq_save = "".join([seq[i] for i in popular_positions])
                print(seq_save)
                f1.write(f"{uniprot}\t{seq_save}\n")
    fin = open('logo.cln')
    seqs = weblogo.read_seq_data(fin)
    logodata = weblogo.LogoData.from_seqs(seqs)
    logooptions = weblogo.LogoOptions()
    logoformat = weblogo.LogoFormat(logodata, logooptions)
    jpeg = weblogo.logo_formatter.jpeg_formatter(logodata, logoformat)
    with open(save_file, "ab") as a:
        a.write(jpeg)

    print(jpeg)


file = ["../data/final_file_aligned.csv",
        "../data/final_results/divided_groups_without_gene_redundancy_fish_Q_aligned.csv",
        f"../data/final_results/divided_groups_without_gene_redundancy_Chondrichtlhyes_aligned.csv",
        "../data/final_results/divided_groups_without_gene_redundancy_Birds_aligned.csv",
        f"../data/final_results/divided_groups_without_gene_redundancy_Amphibians_aligned.csv",
        f"../data/final_results/divided_groups_without_gene_redundancy_reptiles_aligned.csv",
        "../data/final_results/divided_groups_without_gene_redundancy_fish_aligned.csv",
        "../data/final_results/divided_groups_without_gene_redundancy_mammals_aligned.csv",
        "../data/final_results/divided_groups_without_gene_redundancy_Actinopterygii_aligned.csv"]
save_file = ["../data/final_results/logo.jpg",
             "../data/final_results/logo_fish_Q.jpg",
             "../data/final_results/logo_Chondrichtlhyes.jpg",
             "../data/final_results/logo_Birds.jpg",
             "../data/final_results/logo_Amphibians.jpg",
             "../data/final_results/logo_reptiles.jpg",
             "../data/final_results/logo_fish.jpg",
             "../data/final_results/logo_mammals.jpg",
             "../data/final_results/logo_actinopterygii.jpg"]

for e, i in enumerate(file):
    print(i)
    draw_logo(save_file[e], i)
