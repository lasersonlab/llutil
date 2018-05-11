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

from __future__ import division

import re

import pandas as pd
import numpy as np
from scipy.optimize import leastsq


def onebindingsite(x, A, B, C, D):
    """One binding site model

    y = A * x / (B + x) + C * x + D

    A is max
    B is Kd
    C is nonspecific
    D is background
    """
    return A * x / (B + x) + C * x + D



def logistic4(x, A, B, C, D):
    """4-parameter logistic model

    y = D + (A - D) / (1 + ((x / C)**B))

    A is min
    B is Hill coef
    C is inflection
    D is max
    """
    return D + (A - D) / (1 + ((x / C)**B))


def inv_logistic4(y, A, B, C, D):
    return C * ((A - y) / (y - D)) ** (1 / B)


def residual_fn(fn, x, y):
    def residuals(params):
        return (y - fn(x, *params)).ravel()
    return residuals


def load_elisa(path):
    """Loads tab-delim file with ELISA readings

    col `x`: standard concentrations

    col `y1`: standard readings
    col `y2': another standard
    col `yN`: ...

    col `dil`: the dilution multiplier for the sample

    all other columns are particular samples
    """
    df = pd.read_csv(path, sep='\t', header=0)
    y_cols = [c for c in df.columns if re.match('y[0-9]*', c)]
    sample_cols = list(set(df.columns) - set(y_cols) - {'x', 'dil'})
    x = df[['x']].values
    y = df[y_cols].values

    p0 = (0.5,) * 4
    p, cv = leastsq(residual_fn(logistic4, x, y), p0)

    conc = inv_logistic4(df[sample_cols], *p)




import seaborn as sns

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(np.log10(df[:-1]['x']), df[:-1]['y'], 'o')
# ax.plot(np.log10(df['x'].ix[:10]), df['y2'].ix[:10], 'o')
ax.plot(t, z, 'k-')
# ax.plot(np.log10(conc['a']), df[:6]['a'], 'v')
# ax.plot(np.log10(conc['b']), df[:6]['b'], 'v')
# ax.plot(np.log10(conc['c']), df[:6]['c'], 'v')
fig.show()
