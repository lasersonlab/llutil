# Copyright 2017 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from Bio import pairwise2
from Bio import SeqIO



ref = 'AGGCCGCCTTTAGCATGGGATCGATTACGTATAGGATTCGATCGATCGATGAGATCGACG'
s1 = 'CGATGAGATCGA'
s2 = 'GGGGTATAAAAAATGCATTCGAGAGGCCGCCTTTAGCATGGGATCGATTA'


aln1 = pairwise2.align.globalms(ref, s1, 2, -1, -5, -1, penalize_end_gaps=False)[0]
aln2 = pairwise2.align.globalms(ref, s2, 2, -1, -5, -1, penalize_end_gaps=False)[0]

merge_global_alns(aln1[0], [aln1[1]], aln2[0], [aln2[1]])



-----------------------AGGCCGCCTTTAGCATGGGATCGATTACGTATAGGATTCGATCGATCGATGAGATCGACG
---------------------------------------------------------------------CGATGAGATCGA--
GGGGTATAAAAAATGCATTCGAGAGGCCGCCTTTAGCATGGGATCGATTA---------------------------------




from glob import glob

for ab1_file in glob('*.ab1'):
    sr = SeqIO.read(ab1_file, 'abi-trim')



def merge_global_alns(gapped_ref1, gapped_alns1, gapped_ref2, gapped_alns2):
    """Merge pairwise alignments to the same reference

    gapped_ref1 and gapped_alns1 represent the same underlying reference
    sequence to which other sequences are aligned to.

    gapped_alns1 and gapped_alns2 are lists of gapped sequences, each one
    corresponding to an alignment with the respective gapped reference.

    returns merged_ref which is a new gapped version of the reference and
    merged_alns, which is a list containing all the sequences from gapped_alns1
    and gapped_alns1. All sequences are now gapped correspondingly.
    """
    merge_aln = pairwise2.align.globalxx(gapped_ref1, gapped_ref2)[0]
    assert merge_aln[0] == merge_aln[1]
    merged_ref = merge_aln[0]
    merged_alns1 = []
    merged_alns2 = []
    j1 = 0
    j2 = 0
    for i in range(len(merged_ref)):
        if merged_ref[i] == gapped_ref1[j1]:
            merged_alns1.append([aln[j1] for aln in gapped_alns1])
            j1 += 1
        else:
            merged_alns1.append(['-' for aln in gapped_alns1])
        if merged_ref[i] == gapped_ref2[j2]:
            merged_alns2.append([aln[j2] for aln in gapped_alns2])
            j2 += 1
        else:
            merged_alns2.append(['-' for aln in gapped_alns2])
    merged_alns = [''.join(parts) for parts in zip(*merged_alns1)]
    merged_alns += [''.join(parts) for parts in zip(*merged_alns2)]
    return (merged_ref, merged_alns)













def oligo_gen(seq, min_len, max_len):
    """Generate all possible oligos from seq with length constraints

    seq is Bio.Seq.Seq or string
    """
    for i in range(len(seq) - min_len):
        for j in range(min_len, max_len + 1):
            oligo = seq[i:i + j]
            if len(oligo) == j:
                yield oligo


def dna_mutation_gen(seq):
    """Generate all possible point mutations from DNA seq

    seq is Bio.Seq.Seq

    Does not respect case of letters
    """
    letters = seq.alphabet.letters
    for i in range(len(seq)):
        for letter in letters:
            if letter != seq[i].upper():
                yield seq[:i] + letter + seq[i + 1:]


def inosine_gen(seq):
    """Generate all single inosine mutations in seq

    seq is a Bio.Seq.Seq or str

    Does not respect alphabets
    """
    compat = set('GAT')
    for i in range(len(seq)):
        if seq[i].upper() in compat:
            yield seq[:i] + 'I' + seq[i + 1:]


def random_dna_seq(length):
    s = ''.join([random.choice('ACGT') for _ in range(length)])
    return Seq(s, unambiguous_dna)


DNAProp = namedtuple('DNAProp', ['na', 'gc', 'Tm', 'ss_dG', 'self_hyb_dG'])


def dna_properties_batch(seqs):
    """Return a list of namedtuple with some DNA properties

    seqs is a list[Bio.Seq.Seq or str] representing DNA sequences
    """
    seqs = [str(seq) for seq in seqs]
    gcs = [GC(seq) for seq in seqs]
    Tms = melting_temp(seqs)
    ss_dGs = hybrid_ss_min(seqs)
    self_hyb_dGs = [r[0] for r in hybrid_min(seqs, seqs)]
    return [DNAProp(*tup) for tup in zip(seqs, gcs, Tms, ss_dGs, self_hyb_dGs)]


def dna_properties(seq):
    """Return a namedtuple with some DNA properties

    seq is a Bio.Seq.Seq or str representing DNA sequence
    """
    return dna_properties_batch([seq])[0]


ProtProp = namedtuple(
    'ProtProp',
    ['aa', 'gravy', 'aromaticity', 'isoelectric_point', 'instability', 'aa_counts'])


def protein_properties(seq):
    """Return a tuple with some protein biochemical properties

    seq is a Bio.Seq.Seq or str representing protein sequence
    """
    pa = ProteinAnalysis(seq)

    aa_counts = pa.count_amino_acids()
    arom = pa.aromaticity()
    isoelec = pa.isoelectric_point()
    try:
        instability = pa.instability_index()
    except KeyError:
        instability = None
    try:
        gravy = pa.gravy()
    except KeyError:
        gravy = None

    return ProtProp(aa=str(seq),
                    gravy=gravy,
                    aromaticity=arom,
                    isoelectric_point=isoelec,
                    instability=instability,
                    aa_counts=aa_counts)
