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

import re
import os.path as osp
from collections import namedtuple


def to_bytes(bytes_or_str):
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def fastx_stem(path):
    m = re.match('(.+)(?:\.fast[aq]|\.fna|\.f[aq])(?:\.gz)?$', osp.basename(path))
    if m is None:
        raise ValueError(
            'Path {} does not look like a fast[aq] file'.format(path))
    return m.group(1)


def parse_illumina_fastq_name(path):
    """Parse Illumina fastq file name"""
    stem = fastx_stem(path)
    m = re.match('(.*)_S(\d+)_L(\d+)_([RI][12])_001', stem)
    IlluminaFastq = namedtuple(
        'IlluminaFastq', ['sample', 'sample_num', 'lane', 'read', 'path'])
    return IlluminaFastq(sample=m.group(1),
                         sample_num=int(m.group(2)),
                         lane=int(m.group(3)),
                         read=m.group(4),
                         path=path)
