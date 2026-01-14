# Copyright (C) 2026 Zhaoqin Xu, Liu Gan, Mingjie Wei
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import yaml
from pathlib import Path
import basicfunc

# Create an argument parser
parser = argparse.ArgumentParser()

# Add the filename argument
parser.add_argument('filename', type=str, help='The path to the file')

# Parse the command line arguments
args = parser.parse_args()

# Get the filename argument value
filename = args.filename

# Read the file and process the data
with open(filename) as f:
    data = [basicfunc.Dealine(line.strip()) for line in f if line.strip()]

types = {}
flag = None

for line in data:
    if line:
        if type(line[0]) != int:
            flag = line[0]+'s'
            types[flag] = {}
        else:
            types[flag][line[1]] = int(line[0])

# Write the types dictionary to a YAML file
yaml_data = yaml.dump(types, default_flow_style=False)

# Write the YAML data to the file
output_file = Path("testdata/types.yaml")
output_file.write_text(yaml_data)