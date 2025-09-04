# Evolution of Abeta1-42 peptides

## Dataset preparation

To identify proteins containing amyloid beta peptide, we searched for similar proteins using BLAST in UniprotKB and NCBI
non-redundant databases with human abeta1-42. Subsequently, we identified proteins with an analysed fragment containing
histidines. We retrieved the proteins by BLAST and we generated a multiple sequence alignment (MSA) using the MAFFT
method. We localised fragments aligned to ABeta1-42 and filtered out proteins with ABeta6-14 regions shorter than 5
amino acids.
In Actinopterygii, ABeta is encoded by two genes, appa and appb, whereas other organisms have only one gene encoding
ABeta. Therefore, to perform the analysis, we identified one representative ABeta1-42 sequence produced by genes and
excluded isoforms. To find them, we proceeded as follows. We divided the proteins per organism and for each of them we
made an MSA with the human sequence. For each organism, we identified the ABeta1-42 localisations in the MSA and grouped
the same peptides. For each gene, we selected the peptides that were most similar to the human sequence, based on the
Levensthein distance.
For Actinopterygii, we assigned genes appa or appb that encode proteins to each protein. To find the genes, we used the
NCBI, Ensembl and UniprotKB databases. We then executed CD-HIT to group proteins that did not have assigned genes. We
assigned genes to peptides that clustered with peptides with known genes. For instance, a peptide similar to the product
of appa was automatically included in the appa cluster. After thorough examination of parameter combinations, we have
selected a restrictive similarity procedure (similarity: 90%, 95% and 97%). Finally, we selected organisms that had
assigned proteins expressed both genes. If there was more than one peptide variant per organism and gene, we selected
only the most similar one to human ABeta1-42.

We calculated the number of unique proteins containing histidines in MSA at positions 6, 13, and 14 for each animal
group. Additionally, we enumerated the specific mutations at these positions within each group (see Table S1 in
publication). Furthermore, we generated phylogenetic trees and WebLogos to visualize the most significant findings.

## Tree and logo preparation

We performed analyses on the final dataset. A phylogenetic tree was constructed using the ETE3 Python package, based on
the ncbi taxonomy(see Figure 1 in publication. The tree leaves represented animal groups (fish, amphibians, mammals,
reptiles, and birds) and contained counts of ABeta1-42 peptides with histidine positions at residues 6, 13, and 14. For
the group of fish known as Actinopterygii, there are counts for the appa gene, orthologous to the human app gene. The
presence of two leaves with bony and cartilaginous fishes in the group of fishes indicates significant differences in
the sequence peptides. Additionally, the MSA was generated with MAFFT, and the most diverse fragments (residues 1-16)
were visualised with WebLogos for each animal group. WebLogo representations were also provided for the appb gene of
Actinopterygii and the 17-42 peptide fragments were presented for each animal group (see Table S2). Furthermore, the
histidine positions in the peptides from each organism were enumerated and their detailed composition was analysed (see
Table S1).

## Setup and Run code

To prepare the environment for running the analysis, the following must be installed: Python 3.9, CD-HIT and MAFFT.
Additionally,
Installation of the Python 3 packages specified in the requirements.txt file is also required.

First, the region of interest should be searched for in the UniProt and non-redundant databases (the results of which
will be published).
In the data folder. Next, run non-redundant protein selection with the following command:

```
python3 select_nonredundant_proteins.py --fasta_uniprot_file <file with proteins from UniProt in fasta format> --fasta_refseq_file <file with proteins from RefSeq in fasta format> --result_folder <folder that will contain results>
```

Example values for the parameters can be:

- fasta_uniprot_file: "./data/abeta_uniprot.txt"
- fasta_refseq_file: "
  ./data/abeta_refseq.txt"
- the result folder: "./result".

The next commands will calculate the tables or plot the figures that were used in the publication.

```
python3 plot_logo.py --final_file <path> --TaxIds <list> --result_folder <path>
python3 plot_tree.py --final_file <path> --organisms <dict> --aa_positions <dict> --H_positions_file <path>
python3 count_position_z_calue_table.py --final_file <path>  --organisms <lis> --result_folder <path>
python3 plot_Z_hist.py --final_file <path> --files <dict> 
python3 general_mutation_table.py --final_file <path> --aminoacids_positions <dict> --result_file <path>
python3 specific_mutations_table.py --final_file <path> --TaxIds <list> --aminoacids_positions <list> --result_file <path>
```

The final file is generated using the script select_non_redundant_proteins and is named final_file_aligned.csv. This
file is placed in the selected results folder. For the example data, this should be './result/final_file_aligned.csv'.
The commands have parameters with proposed default values that the user can run without specifying.
