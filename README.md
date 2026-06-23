<p align="right">
  <a href="README.md">English</a> | <a href="README.zh-CN.md">中文</a>
</p>

# Welcome to IPS!

<p align="center">
  <img src="assets/ips-logo.svg" alt="IPS logo" width="620">
</p>

# Quick Start

Before getting started, please make sure that you have installed Python 3.5 or later versions and LAMMPS version 3 Mar 2020 or later.

1. Confirm your Python environment

   ```shell
   python --version
   ```

   Confirm whether your Python environment is installed correctly. You should see the following output:

   ```shell
   $ python --version
   Python 3.10.6
   ```

2. Confirm your LAMMPS environment


   ```shell
   lmp -h
   ```

   If there is no error, it means that the LAMMPS environment is normal.

3. Rename the LAMMPS module as `data.lmps` and put it in the project root directory.

4. Put the related input and settings modules into the `./config` folder.

5. Install the required modules by executing the following command in the project root directory:

   ```
   pip install -r requirements.txt
   ```

   Then execute:

   ```shell
   python Simulator.py
   ```

   The program will run.

# Configuration Guide

The configuration file uses YAML format.

Please note that YAML files have strict spacing requirements. For example, the following format is incorrect:

``` yaml
minfile:"config/min.in" # Incorrect format: no space between minfile: and "config/min.in"
```

Correct format:

``` yaml
minfile: "config/min.in" # Correct format: space between minfile: and "config/min.in"
```

Please edit according to YAML format specifications.

## Configuration Parameters

### Basic Program Settings

`minfile`:
- Description: Path to the energy minimization input file
- Value type: String
- Example:
  ```yaml
  minfile: "config/min.in"
  ```

`mdfile`:
- Description: Path to the molecular dynamics input file
- Value type: String
- Example:
  ```yaml
  mdfile: "config/md.in"
  ```

`bondsNum`:
- Description: Current number of formed bonds. Set to `0` to start a new calculation. If set to a value > `0`, the calculation will continue based on data in the `outPath` directory.
- Example:
  ```yaml
  bondsNum: 0
  ```

`bondLimits`:
- Description: Maximum number of bonds that can be formed during the calculation. When this limit is reached, the energy minimization calculation will stop.
- Example:
  ```yaml
  bondLimits: 500
  ```

`tryNum`:
- Description: Number of attempts to find a valid bond pair within the current bonding cycle. If no valid pair is found, the program will try again until a bondable pair is found or the maximum number of attempts is reached.
- Example:
  ```yaml
  tryNum: 1
  ```

`trylimit`:
- Description: Maximum number of attempts per bonding cycle to find a valid bond pair. If this limit is exceeded without success, the entire calculation will terminate.
- Example:
  ```yaml
  trylimit: 140
  ```

`lammps`:
- Description: Command to execute LAMMPS from the command line.
- Example:
  ```yaml
  lammps: "mpirun -np 32 lmp_oneapi"
  ```

`outPath`:
- Description: Directory path for storing calculation results from each iteration.
- Example:
  ```yaml
  outPath: "outdirs"
  ```

`lmprundir`:
- Description: Directory path where LAMMPS runs.
- Example:
  ```yaml
  lmprundir: "lammps"
  ```

### Bonding Condition Settings

These settings define the constraints for bond formation.

`link`:
- Description: Defines the atom types to be bonded. Each entry consists of pairs of atom types. The atom type codes must be defined in the molecular system settings. Supports multiple link groups for different bonding scenarios.
- Example:
  ```yaml
  link:
    -
      - nh,n
      - c,c_1
    -
      - tnv,tns
      - c,c_1
  ```

`deletelink`:
- Description: Atom types to be deleted when forming bonds (e.g., leaving groups).
- Example:
  ```yaml
  deletelink:
    - hn
    - cl
    - thn
  ```

`cutoff`:
- Description: Maximum distance between atoms that can form a bond, in Angstroms (Å).
- Example:
  ```yaml
  cutoff: 3.5
  ```

`charge`:
- Description: Virtual charges to be assigned to the two bonding atoms. The first and second values correspond to the atom types defined in the same positions in `link`.
- Example:
  ```yaml
  charge:
    - 0
    - 0
  ```

`intra`:
- Description: The first parameter defines whether to check if the two potentially bondable atoms are adjacent atoms on the same molecule. The second parameter defines the maximum number of atoms between them for this adjacency relationship.
- Example:
  ```yaml
  intra:
    - true
    - 1
  ```

### Equilibrium Settings

These settings are used when equilibrium calculations are enabled.

`equilibrium`:
- Description: Whether to perform equilibrium calculations during the simulation.
- Value type: Boolean
- Example:
  ```yaml
  equilibrium: true
  ```

When `equilibrium: true`, the following parameters are required:

`upperPlate`:
- Description: Atom type name for the upper plate.
- Value type: String
- Example:
  ```yaml
  upperPlate: "up"
  ```

`upperSolvent`:
- Description: List of atom type numbers for the upper solvent.
- Value type: List of integers
- Example:
  ```yaml
  upperSolvent:
    - 11
    - 12
  ```

`upperMonomer`:
- Description: LAMMPS data file path for the upper monomer.
- Value type: String
- Example:
  ```yaml
  upperMonomer: "config/tmc.lmps"
  ```

`upperFactor`:
- Description: Ratio factor for comparing upper monomer count with solvent count. The monomer count is multiplied by this factor before comparison.
- Value type: Integer
- Example:
  ```yaml
  upperFactor: 97
  ```

`lowerPlate`:
- Description: Atom type name for the lower plate.
- Value type: String
- Example:
  ```yaml
  lowerPlate: "be"
  ```

`lowerSolvent`:
- Description: List of atom type numbers for the lower solvent.
- Value type: List of integers
- Example:
  ```yaml
  lowerSolvent:
    - 13
    - 14
  ```

`lowerMonomer`:
- Description: LAMMPS data file path for the lower monomer.
- Value type: String
- Example:
  ```yaml
  lowerMonomer: "config/mpd.lmps"
  ```

`lowerFactor`:
- Description: Ratio factor for comparing lower monomer count with solvent count. The monomer count is multiplied by this factor before comparison.
- Value type: Integer
- Example:
  ```yaml
  lowerFactor: 258
  ```

### Atom Type Information

These parameters define the type indices corresponding to `atoms`, `bonds`, `angles`, `dihedrals`, and `impropers` in `data.lmps`. For impropers, the first atom is the central atom.

Example:

```yaml
atoms:
  c: 1
  c3: 2
  ca: 3
  ha: 4

bonds:
  c3,c3: 1
  c3,hc: 2
  c3,n: 3
  c3,n3: 4

angles:
  c3,c3,hc: 1
  c3,c3,n3: 2
  c3,c3,n: 3
  hc,c3,hc: 4
  hc,c3,n3: 5

dihedrals:
  ca,ca,c,cl: 1
  ca,ca,c_1,n: 2
  o,c_1,n,c3: 3
  ca,c_1,n,c3: 4
  n3,c3,c3,n: 5

impropers:
  ca,o,c_1,n: 1
  ca,o,c,cl: 2
  ca,ca,ca,ha: 3
```

# Citation

If you use this software in your work, please cite:

```bibtex
@article{ips_2026,
  title = {{IPS}: {An} interfacial polymerization simulator for modeling the highly-cross-linked polymer membranes},
  author = {Liu, Gan and Xu, Zhaoqin and Wei, Mingjie and Huang, Jun and Du, Yan and Wang, Yong},
  journal = {Advanced Membranes},
  volume = {9},
  pages = {100252},
  year = {2026},
  doi = {10.1016/j.advmem.2026.100252}
}
```


# License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See [LICENSE](LICENSE) for details.
