class Protein:
    def __init__(self, uniprot_id, seq, organism_id=None, gene=None, emsembl=None, orf=None):
        self.seq = seq
        self.uniprot_id = uniprot_id
        self.organism_id = organism_id
        # self.right_place_histidine = self.find_histidine()
        self.gene = gene
        self.type = "general"
        self.sequence_redundancy = []
        self.full_sequence = ""
        self.emsembl = emsembl
        self.orf = orf
        self.header = ""
        self.line = ""

    def find_histidine(self, positions):
        is_pos = 0
        for aminoacid, position in positions.items():
            if self.seq[position] == aminoacid:
                is_pos = 1
        if is_pos:
            return True
        else:
            return False

    def get_organism(self, ncbi, typ="general", selected=[7777, 40674, 1489341, 32561, 7894, 8292, 7898, 8782]):
        if typ == "general":
            return self.get_organism_general(ncbi, selected=selected)
        elif typ == "diff":
            return self.get_organism_prec(ncbi)


        elif self.type == "general":
            return self.get_organism_general(ncbi, selected=selected)
        else:
            return self.get_organism_prec(ncbi)

    def get_organism_general(self, ncbi, selected):
        organisms = None
        if self.organism_id != "check":
            org = ncbi.get_lineage(self.organism_id)
            for i in selected:
                if i in org:
                    organisms = i
        else:
            return self.organism_id

        return organisms

    def get_organism_prec(self, ncbi):
        organisms = None
        groups = [7777, 7894, 7899, 7914, 41712, 41705, 123369, 186626, 29140, 32446]
        if self.organism_id != "check":
            org2 = {j: i for j, i in ncbi.get_rank(ncbi.get_lineage(self.organism_id)).items()}
            for i in groups:
                if i in org2:
                    return organisms
            else:
                return self.organism_id
        else:
            return self.organism_id

    def find_by_positions_histidine(self, positions):
        is_pos = {f"{i}{j + 1}": 0 for i, j in positions.items()}
        for aminoacid, position in positions.items():
            if self.seq[position] == aminoacid[0]:
                is_pos[f"{aminoacid}{position + 1}"] = 1
        return is_pos
