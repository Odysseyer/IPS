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

import numpy as np


class atom(object):
    index: int = 0
    type: int = 0
    mol: int = 0
    x: float = 0
    y: float = 0
    z: float = 0
    abbreviation = ""

    def __init__(self, var: list) -> None:
        self.index: int = int(var[0])
        self.mol: int = int(var[1])
        self.type: int = int(var[2])
        self.charge: float = float(var[3])
        self.x: float = float(var[4])
        self.y: float = float(var[5])
        self.z: float = float(var[6])
        self.bonds: list[multi] = []
        self.bonded: list[atom] = []
        self.angles: list[multi] = []
        self.diheds: list[multi] = []
        self.improps: list[multi] = []

    # Get atom type string
    def findtypestr(self, typestruct: dict) -> None:
        for k, v in typestruct.items():
            if v == self.type:
                self.abbreviation: str = k
                return
        raise KeyError(
            f"Not find {self.type} in atom_types , index:{self.index} Mol:{self.mol}"
        )

    @property
    def position(self):
        return np.array([self.x, self.y, self.z])

    def addbonds(self, var: list):
        self.bonds.append(var)

    def __lt__(self, other):
        return self.index < other.index

    def __str__(self) -> str:
        return f"{self.index} {self.abbreviation}"

    def __repr__(self) -> str:
        return f"{self.index} {self.abbreviation}"

    def __eq__(self, value: object) -> bool:
        if self.index==value.index and self.type==value.type and self.mol==value.mol and self.charge==value.charge and self.x==value.x and self.y==value.y and self.z==value.z:
            return True
        else:
            return False


class multi(object):
    def __init__(self, var: list) -> None:
        self.index: int = int(var[0])
        self.type: int = int(var[1])
        self.atoms: list[atom] = var[2:]
        self.atomsI: list[atom] = []

    def LoadAtoms(self, Atoms: list[atom]) -> None:
        for atomindex in self.atoms:
            self.atomsI.append(Atoms[atomindex - 1])

    def update(self) -> None:
        self.atoms = [atom.index for atom in self.atomsI]

    def __lt__(self, other):
        return self.index < other.index
    
    def __eq__(self, value: object) -> bool:
        if self.index==value.index and self.type==value.type and self.atoms==value.atoms:
            return True
        else:
            return False
        
    def __hash__(self) -> int:
        return hash(self.index)

    @property
    def atoms_num(self) -> list[int]:
        return [atom.index for atom in self.atoms]

    @property
    def abbreviation(self) -> str:
        temp = []
        for atom in self.atoms:
            temp.append(atom.abbreviation)
        return ",".join(temp)

    def __str__(self) -> str:
        return f"{self.index} {self.abbreviation}"

    def __repr__(self) -> str:
        return f"{self.index} {self.abbreviation}"

    def __eq__(self, value: object) -> bool:
        if (
            self.index == value.index
            and self.type == value.type
            and self.atoms == value.atoms
        ):
            return True
        else:
            return False


class bond(multi):
    def __init__(self, var) -> None:
        super().__init__(var)


class angle(multi):
    def __init__(self, var: list) -> None:
        super().__init__(var)


class dihedral(multi):
    def __init__(self, var: list) -> None:
        super().__init__(var)


class improper(multi):
    def __init__(self, var: list) -> None:
        super().__init__(var)


class Molecular(object):

    def __init__(self) -> None:
        self.atoms: dict[int, atom] = {}
        self.bonds: dict[int, bond] = {}
        self.angles: dict[int, angle] = {}
        self.dihedrals: dict[int, dihedral] = {}
        self.impropers: dict[int, improper] = {}

    def readMononerfile(self, filename) -> None:
        atoms = {}
        bonds = {}
        angles = {}
        dihedrals = {}
        impropers = {}

        with open(filename, "r") as file:
            lines = file.readlines()

        current_section = None
        for line in lines:
            line = line.strip()
            if line:
                if line.startswith("Atoms"):
                    current_section = "Atoms"
                elif line.startswith("Bonds"):
                    current_section = "Bonds"
                elif line.startswith("Angles"):
                    current_section = "Angles"
                elif line.startswith("Dihedrals"):
                    current_section = "Dihedrals"
                elif line.startswith("Impropers"):
                    current_section = "Impropers"
                elif line.startswith("#") or len(line) == 0:
                    continue
                else:
                    if current_section == "Atoms":
                        atoms[int(line.split()[0])] = atom(
                            list(map(float, line.split()))
                        )
                    elif current_section == "Bonds":
                        temp = list(map(int, line.split()))
                        bonds[temp[0]] = bond(temp[:2] + [atoms[a] for a in temp[2:]])
                    elif current_section == "Angles":
                        temp = list(map(int, line.split()))
                        angles[temp[0]] = angle(temp[:2] + [atoms[a] for a in temp[2:]])
                    elif current_section == "Dihedrals":
                        temp = list(map(int, line.split()))
                        dihedrals[temp[0]] = dihedral(
                            temp[:2] + [atoms[a] for a in temp[2:]]
                        )
                    elif current_section == "Impropers":
                        temp = list(map(int, line.split()))
                        impropers[temp[0]] = improper(
                            temp[:2] + [atoms[a] for a in temp[2:]]
                        )

        self.atoms: dict[int, atom] = atoms
        self.bonds: dict[int, bond] = bonds
        self.angles: dict[int, angle] = angles
        self.dihedrals: dict[int, dihedral] = dihedrals
        self.impropers: dict[int, improper] = impropers

    def update(
        self,
        atom_num: int,
        bond_num: int,
        angle_num: int,
        dihedral_num: int,
        improper_num: int,
        Mols_num: int
    ) -> None:
        new_atoms = {}
        new_bonds = {}
        new_angles = {}
        new_dihedrals = {}
        new_impropers = {}
        
        for key, atom in self.atoms.items():
            atom.index += atom_num
            atom.mol = Mols_num + 1
            new_atoms[atom.index] = atom


        for key, bond in self.bonds.items():
            bond.index += bond_num
            new_bonds[bond.index] = bond

        for key, angle in self.angles.items():
            angle.index += angle_num
            new_angles[angle.index] = angle

        for key, dihedral in self.dihedrals.items():
            dihedral.index += dihedral_num
            new_dihedrals[dihedral.index] = dihedral

        for key, improper in self.impropers.items():
            improper.index += improper_num
            new_impropers[improper.index] = improper

        self.atoms = new_atoms             
        self.bonds = new_bonds
        self.angles = new_angles
        self.dihedrals = new_dihedrals
        self.impropers = new_impropers
       

    @property
    def type(self) -> set:
        return set([atom.type for atom in self.atoms.values()])


