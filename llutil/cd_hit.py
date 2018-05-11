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


from click import group, option, File
import matplotlib as mpl
# preemptively load a matplotlib backend to get around issue with framework
mpl.use('agg')
import skbio

@group(name='cd-hit')
def cd_hit_cli():
    """cd-hit -- commands for working with CD-HIT output"""
    pass


@cd_hit_cli.command()
@option('--representatives', '-r', type=File('r'),
        help='The .fastx output with the cluster representatives')
@option('--clusters', '-c', type=File('r'),
        help='The .clstr file with the cluster definitions')
@option('--output', '-o', type=File('w'), help='Output file')
def collapse(representatives, clusters, output):
    """collapse each cluster into a (size, representative) tuple"""
    # first accumulate the cluster sizes
    cluster_sizes = {}
    first_line = False
    cluster_size = None
    cluster_rep_id = None
    for line in clusters:
        if line.startswith('>'):
            if cluster_rep_id is not None:
                cluster_sizes[cluster_rep_id] = cluster_size
            first_line = True
        elif first_line:
            raw_fields = line.split()
            # some predicates on conditions to test for
            non_asterisk = raw_fields[-1].strip() != '*'
            non_ellipsis = raw_fields[2][-3:] != '...'
            non_gt = raw_fields[2][0] != '>'
            if non_asterisk or non_ellipsis or non_gt:
                raise ValueError(
                    'Unexpected value in first cluster line:\n{}'.format(line))
            cluster_rep_id = raw_fields[2][1:-3]
            cluster_size = 1
            first_line = False
        else:
            cluster_size += 1
    # handle the last line
    cluster_sizes[cluster_rep_id] = cluster_size
    # then pull the associated representative sequences and write out
    format = skbio.io.sniff(representatives)[0]
    kwargs = {} if format == 'fasta' else {'variant': 'illumina1.8'}
    for seq in skbio.io.read(representatives, format=format, **kwargs):
        id_ = seq.metadata['id']
        print('{}\t{}\t{}'.format(cluster_sizes[id_], id_, str(seq)),
              file=output)
