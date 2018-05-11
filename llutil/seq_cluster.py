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

import numpy as np
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage
from jellyfish import damerau_levenshtein_distance


# scipy.spatial.distance.pdist cannot handle strings
def pdist(X, metric):
    m = len(X)
    dm = np.zeros(m * (m - 1) // 2, dtype=np.double)
    k = 0
    for i in range(0, m - 1):
        for j in range(i + 1, m):
            dm[k] = metric(X[i], X[j])
            k += 1
    return dm


seqs = list(SeqIO.parse('/Users/laserson/lasersonlab/airr_primers/isotype-distinguishing/IGHC-CH1.fasta', 'fasta'))
dist_mat = pdist([str(s.seq) for s in seqs], damerau_levenshtein_distance)
clustering = linkage(dist_mat, method='average')
