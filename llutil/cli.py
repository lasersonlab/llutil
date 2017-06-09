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


from click import group

from llutil.cd_hit import cd_hit_cli
from llutil.indexes import make_indexes


@group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """llutil -- Laserson Lab Utilities"""
    pass


cli.add_command(cd_hit_cli)
cli.add_command(make_indexes)
