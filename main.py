COMPLEMENT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
}

GEN_CODE = {
    'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', 'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
    'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V', 'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
    'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A', 'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
    'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
    'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D', 'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
    'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'TGT': 'C', 'CGT': 'R',
    'AGT': 'S', 'GGT': 'G', 'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 'CGA': 'R', 'AGA': 'R',
    'GGA': 'G', 'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

STOP = ['TAG', 'TGA', 'TAA']
START = ['ATG']


def rev_comp(seq):
    rev_seq = ""
    for base in seq:
        rev_seq = COMPLEMENT[base] + rev_seq
    return rev_seq


def orfs(seq):
    o = []

    for i in range(0, len(seq)):
        codon = seq[i:i + 3]
        if codon in START:

            current_orf = []
            j = i
            found_stop = False

            while not found_stop and j < len(seq):
                codon_orf = seq[j:j + 3]
                if codon_orf in STOP:
                    found_stop = True
                else:
                    current_orf.append(codon_orf)

                j += 3

            if found_stop:
                o.append(current_orf)

    return o


if __name__ == '__main__':
    dna = open('seq.txt', 'r').read().replace("\n", "")
    dna_rev = rev_comp(dna)

    seqs = [dna, dna_rev]

    # print(seqs)

    orfs_codons = []
    for seq in seqs:
        orfs_codons += orfs(seq)

    # print(orfs_codons)

    peps = [''.join(map(lambda codon: GEN_CODE[codon], orf)) for orf in orfs_codons]

    # print(peps)

    unique_peps = set()
    for pep in peps:
        unique_peps.add(pep)

    print("===== DNA input =====")
    print(dna)
    print()
    print("======= ORFS ========")
    with open('out.txt', 'w') as f:
        for upep in unique_peps:
            print(upep)
            f.write(upep + '\n')
