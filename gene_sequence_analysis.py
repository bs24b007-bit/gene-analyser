#!/usr/bin/env python3
"""Gene Sequence Analysis Tool

Usage:
  python gene_sequence_analysis.py --sequence ATGCGT...
  python gene_sequence_analysis.py --file sequence.fasta
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

CODON_TABLE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

NUCLEOTIDE_WEIGHTS: Dict[str, float] = {
    'A': 313.21,
    'C': 289.18,
    'G': 329.21,
    'T': 304.2,
}

AMINO_ACID_WEIGHTS: Dict[str, float] = {
    'A': 89.09, 'R': 174.2, 'N': 132.12, 'D': 133.1,
    'C': 121.15, 'Q': 146.15, 'E': 147.13, 'G': 75.07,
    'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
    'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
    'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
    '*': 0.0,
}

COMPLEMENT: Dict[str, str] = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'N': 'N',
}

START_CODON = 'ATG'
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def normalize_sequence(sequence: str) -> str:
    return re.sub(r'[^ACGTacgt]', '', sequence.upper())


def parse_fasta(file_path: str) -> Tuple[str, str]:
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Path not found: {file_path}")

    header = ''
    sequence_lines: List[str] = []

    with open(file_path, 'r', encoding='utf-8') as fasta:
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if not header:
                    header = line[1:].strip()
                continue
            sequence_lines.append(line)

    if not sequence_lines:
        raise ValueError('No sequence found in FASTA file.')

    return header or os.path.basename(file_path), normalize_sequence(''.join(sequence_lines))


def gc_content(sequence: str) -> float:
    sequence = normalize_sequence(sequence)
    if not sequence:
        return 0.0
    gc_count = sum(1 for base in sequence if base in ('G', 'C'))
    return 100.0 * gc_count / len(sequence)


def reverse_complement(sequence: str) -> str:
    sequence = normalize_sequence(sequence)
    return ''.join(COMPLEMENT.get(base, 'N') for base in reversed(sequence))


def translate_dna(sequence: str, frame: int = 1) -> str:
    seq = normalize_sequence(sequence[frame - 1:])
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)


def find_orfs(sequence: str, min_length: int = 90) -> List[Dict[str, object]]:
    sequence = normalize_sequence(sequence)
    orfs: List[Dict[str, object]] = []
    strands = [('+', sequence), ('-', reverse_complement(sequence))]

    for strand_label, strand_seq in strands:
        for frame in range(3):
            i = frame
            while i <= len(strand_seq) - 3:
                codon = strand_seq[i:i + 3]
                if codon == START_CODON:
                    j = i
                    while j <= len(strand_seq) - 3:
                        current = strand_seq[j:j + 3]
                        if current in STOP_CODONS:
                            orf_sequence = strand_seq[i:j + 3]
                            if len(orf_sequence) >= min_length:
                                orfs.append({
                                    'strand': strand_label,
                                    'frame': frame + 1,
                                    'start': i + 1,
                                    'end': j + 3,
                                    'length': len(orf_sequence),
                                    'dna': orf_sequence,
                                    'protein': translate_dna(orf_sequence),
                                })
                            break
                        j += 3
                i += 3

    orfs.sort(key=lambda x: x['length'], reverse=True)
    return orfs


def molecular_weight_dna(sequence: str) -> float:
    sequence = normalize_sequence(sequence)
    return sum(NUCLEOTIDE_WEIGHTS.get(base, 0.0) for base in sequence)


def molecular_weight_protein(sequence: str) -> float:
    return sum(AMINO_ACID_WEIGHTS.get(residue, 0.0) for residue in sequence)


def format_orf(orf: Dict[str, object]) -> str:
    return (
        f"Strand {orf['strand']} frame {orf['frame']}: "
        f"{orf['start']}..{orf['end']} ({orf['length']} bp)\n"
        f"  DNA: {orf['dna']}\n"
        f"  Protein: {orf['protein']}"
    )


def summarize(sequence: str, name: Optional[str] = None) -> str:
    seq = normalize_sequence(sequence)
    if not seq:
        return 'No valid DNA sequence provided.'

    rc = reverse_complement(seq)
    gc = gc_content(seq)
    dna_weight = molecular_weight_dna(seq)
    protein = translate_dna(seq)
    protein_weight = molecular_weight_protein(protein)
    orfs = find_orfs(seq)

    lines = [
        f"Sequence name: {name or 'input'}",
        f"Length: {len(seq)} bp",
        f"GC content: {gc:.2f}%",
        f"Reverse complement: {rc}",
        f"DNA molecular weight: {dna_weight:.2f} g/mol",
        f"Protein translation (frame 1): {protein}",
        f"Protein molecular weight: {protein_weight:.2f} Da",
        f"Open reading frames found: {len(orfs)}",
    ]

    for index, orf in enumerate(orfs[:10], start=1):
        lines.append(f"\nORF {index}:")
        lines.append(format_orf(orf))

    return '\n'.join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description='Gene Sequence Analysis Tool')
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--sequence', '-s', help='Enter a DNA sequence directly')
    group.add_argument('--file', '-f', help='Path to a FASTA file')
    parser.add_argument('--min-orf-length', type=int, default=90, help='Minimum ORF length in bp')

    args = parser.parse_args()

    sequence = ''
    name = None

    if args.file:
        name, sequence = parse_fasta(args.file)
    elif args.sequence:
        sequence = args.sequence.strip()
        name = 'input sequence'
    else:
        print('Enter DNA sequence or FASTA file path. Press Enter when done.')
        user_input = sys.stdin.read().strip()
        if os.path.isfile(user_input):
            name, sequence = parse_fasta(user_input)
        else:
            sequence = user_input
            name = 'input sequence'

    if not sequence:
        print('No sequence provided. Use --sequence or --file, or pipe a sequence into the tool.')
        sys.exit(1)

    sequence = normalize_sequence(sequence)
    output = summarize(sequence, name=name)
    print(output)


if __name__ == '__main__':
    main()
