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


from subprocess import check_call
import sys

from click import command, option, argument, Path

from llutil.util import fastx_stem


@command(name='make-indexes')
@argument('fasta', type=Path(exists=True, dir_okay=False))
@option('--name', '-n',
        help=('A name for the index. Uses the stem of the fasta file if '
              'unspecified.'))
@option('--output-dir', '-o', type=Path(exists=True, file_okay=False),
        help=('Output directory to write indexes. Uses the current working '
              'directory if unspecified.'))
@option('--dry-run', '-d', is_flag=True,
        help='Dry run; print commands that will be run.')
def make_indexes(fasta, name, output_dir, dry_run):
    """make-indexes -- create indexes from fasta file

    Does not accept gzipped input.

    Constructs indexes for bowtie, bowtie2, bwa, kallisto, all of which must be
    available on the PATH.
    """
    if name is None:
        name = fastx_stem(fasta)
    if output_dir is None:
        output_dir = '.'
    if fasta.endswith('.gz'):
        print('Does not support gzipped fasta files; please gunzip first.')
        sys.exit()

    bowtie_cmd = (
        'mkdir {output_dir}/bowtie && bowtie-build -q {fasta} {output_dir}/bowtie/{name}').format(
            output_dir=output_dir, name=name, fasta=fasta)
    bowtie2_cmd = (
        'mkdir {output_dir}/bowtie2 && bowtie2-build -q {fasta} {output_dir}/bowtie2/{name}').format(
            output_dir=output_dir, name=name, fasta=fasta)
    bwa_cmd = (
        'mkdir {output_dir}/bwa && bwa index -p {output_dir}/bwa/{name} {fasta}').format(
            output_dir=output_dir, name=name, fasta=fasta)
    kallisto_cmd = (
        'mkdir {output_dir}/kallisto && kallisto index -i {output_dir}/kallisto/{name}.idx {fasta}').format(
            output_dir=output_dir, name=name, fasta=fasta)

    if dry_run:
        print(bowtie_cmd)
        print(bowtie2_cmd)
        print(bwa_cmd)
        print(kallisto_cmd)
    else:
        check_call(bowtie_cmd, shell=True)
        check_call(bowtie2_cmd, shell=True)
        check_call(bwa_cmd, shell=True)
        check_call(kallisto_cmd, shell=True)
