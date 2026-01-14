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

class instruction:
    link = []
    charge = []
    cutoff = []
    bond = []
    align = []
    alignConfig = []

    def __init__(self, data: dict) -> None:
        self.link = data.get("link", None)
        self.link[0] = list(map(lambda x: x.split(","), self.link[0]))
        if self.link.__len__()==2:
            self.link[1] = list(map(lambda x: x.split(","), self.link[1]))
        self.deletelink = data.get("deletelink", None)
        self.charge = data.get("charge", None)
        self.cutoff = data.get("cutoff", None)
        self.bond = data.get("bond", None)
        self.intra = data.get("intra", None)
        self.align = data.get("align", None)
        self.alignConfig = data.get("alignConfig", None)

        self.equilibrium: bool = data.get("equilibrium", False)  # Whether to perform equilibrium calculation

        if self.equilibrium:
            self.upperPlate: str = data.get(
                "upperPlate", None
            )  # Upper plate, name of atom types composing upper plate
            self.upperSolvent: list[int] = data.get(
                "upperSolvent", None
            )  # Upper plate solvent, list of atom type numbers for upper solvent
            self.upperMonomer: str = data.get(
                "upperMonomer", None
            )  # Upper plate monomer, path to upper monomer's LAMMPS data file
            self.upperFactor: int = data.get(
                "upperFactor", 0
            )  # Ratio of upper monomer to upper solvent, multiply monomer count by this ratio and compare with solvent count

            self.lowerPlate: str = data.get(
                "lowerPlate", None
            )  # Lower plate, name of atom types composing lower plate
            self.lowerSolvent: list[int] = data.get(
                "lowerSolvent", None
            )  # Lower plate solvent, list of atom type numbers for lower solvent
            self.lowerMonomer: str = data.get(
                "lowerMonomer", None
            )  # Lower plate monomer, path to lower monomer's LAMMPS data file
            self.lowerFactor: int = data.get(
                "lowerFactor", 0
            )  # Ratio of lower monomer to lower solvent, multiply monomer count by this ratio and compare with solvent count


class types:
    atom_Type: dict = {}
    bond_Type: dict = {}
    angle_Type: dict = {}
    dihedral_Type: dict = {}
    improper_Type: dict = {}

    def __init__(self, data) -> None:
        # Swap key-value pairs
        data["atoms"].update({value: key for key, value in data["atoms"].items()})
        data["bonds"].update({value: key for key, value in data["bonds"].items()})
        data["angles"].update({value: key for key, value in data["angles"].items()})
        data["dihedrals"].update(
            {value: key for key, value in data["dihedrals"].items()}
        )
        data["impropers"].update(
            {value: key for key, value in data["impropers"].items()}
        )
        self.atom_Type: dict = data["atoms"]
        self.bond_Type: dict = data["bonds"]
        self.angle_Type: dict = data["angles"]
        self.dihedral_Type: dict = data["dihedrals"]
        self.improper_Type: dict = data["impropers"]


