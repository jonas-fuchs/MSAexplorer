AMBIG_NUCS = {
        'DNA': {
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'G', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
        },
        'RNA': {
        'R': ['A', 'G'],
        'Y': ['C', 'U'],
        'S': ['G', 'C'],
        'W': ['A', 'U'],
        'K': ['G', 'U'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'U'],
        'D': ['A', 'G', 'U'],
        'H': ['A', 'G', 'U'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'U']
        }
    }

complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',  # R = A/G, Y = C/T, S = G/C, W = A/T
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',  # K = G/T, M = A/C, B = C/G/T, V = A/C/G
        'D': 'H', 'H': 'D', 'N': 'N', '-': '-',  # D = A/G/T, H = A/C/T, N = any base, '-' = gap
        'U': 'A'                                 # RNA
    }


amino_acids = ['A', 'R', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

start_codons = {
        'DNA': ['ATG'],
        'RNA': ['AUG'],
}

stop_codons = {
        'DNA': ['TAG', 'TGA', 'TAA'],
        'RNA': ['UAG', 'UGA', 'UAA'],
}
