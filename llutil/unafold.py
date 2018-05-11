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

import math
from subprocess import run, PIPE

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna


"""UNAFold wrappers/utils

Download UNAFold from here:
http://homepages.rpi.edu/~zukerm/download/UNAFold_download.html
"""


def hybrid_ss_min(seqs, na='RNA', tmin=37, tinc=1, tmax=37, sodium=1,
                  magnesium=0):
    """Run hybrid-ss-min on a list of sequences

    seqs is a list of Bio.Seq.Seq objects or strings

    returns [deltaG]
    """
    cmd = ('hybrid-ss-min --stream --NA {} --tmin {} --tinc {} --tmax {} '
           '--sodium {} --magnesium {}').format(
               na, tmin, tinc, tmax, sodium, magnesium)
    input_string = ''.join(['{};\n'.format(str(s)) for s in seqs])
    cp = run(cmd, shell=True, input=input_string, stdout=PIPE,
             universal_newlines=True, check=True)
    results = [float(line) for line in cp.stdout.strip().split('\n')]
    if len(results) != len(seqs):
        raise ValueError('num input seqs does not match num of results')
    return results


def hybrid_min(seqs1, seqs2, na='RNA', tmin=37, tinc=1, tmax=37, sodium=1,
               magnesium=0):
    """Run hybrid-min on pairs of sequences

    seqs1 and seqs2 are both lists of Bio.Seq.Seq objects or strings

    returns [(deltaG, deltaH, deltaS)]
    """
    if len(seqs1) != len(seqs2):
        raise ValueError('seqs1 and seqs2 must be the same length')
    cmd = ('hybrid-min --stream --NA {} --tmin {} --tinc {} --tmax {} '
           '--sodium {} --magnesium {}').format(
               na, tmin, tinc, tmax, sodium, magnesium)
    input_string = ''.join(['{};\n{};\n'.format(str(s1), str(s2))
                            for s1, s2 in zip(seqs1, seqs2)])
    cp = run(cmd, shell=True, input=input_string, stdout=PIPE,
             universal_newlines=True, check=True)
    process_line = lambda line: tuple([float(f) for f in line.split('\t')])
    results = [process_line(line) for line in cp.stdout.strip().split('\n')]
    if len(results) != len(seqs1):
        raise ValueError('num input seqs does not match num of results')
    return results


def melting_temp(seqs, oligo_conc=0.25, sodium=0.05, magnesium=0):
    """Get Tm of list of seqs using UNAFold calculation

    seqs is a list of Bio.Seq.Seq or strings
    oligo_conc is in uM
    ion conc is in M

    For Q5, use sodium=0.1, magnesium=0.002

    returns [Tm]
    """
    revcomp = lambda s: str(Seq(str(s), unambiguous_dna).reverse_complement())
    rcs = [revcomp(s) for s in seqs]
    energies = hybrid_min(seqs, rcs, na='DNA', sodium=sodium, magnesium=magnesium)
    # Annu. Rev. Biophys. Biomol. Struct. 2004. 33:415–40
    # doi: 10.1146/annurev.biophys.32.110601.141800
    R = 1.9872036
    molar_conc = oligo_conc * 1e-6
    tm = lambda dH, dS: dH * 1000 / (dS + R * math.log(molar_conc)) - 273.15
    return [tm(dH, dS) for dG, dH, dS in energies]
