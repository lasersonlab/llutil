# Copyright 2016 Uri Laserson
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

import random

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna


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
