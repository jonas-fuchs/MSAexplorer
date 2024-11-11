POSSIBLE_CHARS = [
        'A', 'R', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*', 'X',  #AS
        'A', 'T','U', 'C', 'G', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V','N',  #RNA/DNA
        '-'  #GAP
]

AMBIG_CHARS = {
        'DNA': {
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
        '-': ['A', 'C', 'G', 'T'],
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
        'H': ['A', 'C', 'U'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'U'],
        '-': ['A', 'C', 'G', 'T'],
        },
        'AS': {
        'X': ['A', 'R', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*'],
        '-': ['A', 'R', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        }
    }

COMPLEMENT = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',  # R = A/G, Y = C/T, S = G/C, W = A/T
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',  # K = G/T, M = A/C, B = C/G/T, V = A/C/G
        'D': 'H', 'H': 'D', 'N': 'N', '-': '-',  # D = A/G/T, H = A/C/T, N = any base, '-' = gap
        'U': 'A'                                 # RNA
    }

START_CODONS = {
        'DNA': ['ATG'],
        'RNA': ['AUG'],
}

STOP_CODONS = {
        'DNA': ['TAG', 'TGA', 'TAA'],
        'RNA': ['UAG', 'UGA', 'UAA'],
}
