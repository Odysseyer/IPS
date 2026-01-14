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

import basicfunc
import Polymer

class readinfile:
    def __init__(self, filename) -> None:
        with open(filename, "r", encoding="utf-8") as f:
            indata = f.readlines()
        tags = ["read_data", "displace_atoms", "group", "variable", "change_box"]
        data = {tag: [] for tag in tags}
        for line in indata:
            # If the beginning of this line matches one of our tags, store it
            tag, parts = self.parse_line(line)
            if tag in tags:
                data[tag].append(parts)

    def parse_line(line: str):
        parts = line.split()
        return parts[0], parts[1:]


class CoeffsInfoclass:  # Used to store Coeffs information
    def __init__(self) -> None:
        self.Angle_Coeffs: list = []
        self.Bond_Coeffs: list = []
        self.Dihedral_Coeffs: list = []
        self.Masses: list = []
        self.Pair_Coeffs: list = []
        self.Improper_Coeffs: list = []


class MolsysInfoclass:  # Used to store system information
    angle_type: int = 0
    angles_num: int = 0
    atom_type: int = 0
    atoms_num: int = 0
    bond_type: int = 0
    bonds_num: int = 0
    dihedral_type: int = 0
    dihedrals_num: int = 0
    improper_type: int = 0
    impropers_num: int = 0
    xhi: float = 0
    xlo: float = 0
    yhi: float = 0
    ylo: float = 0
    zhi: float = 0
    zlo: float = 0


class lammps:  # LAMMPS class
    def __init__(self, filename) -> None:
        temp = readlammps(filename)
        self.Atoms: dict[int, Polymer.atom] = {}
        self.Bonds: dict[int, Polymer.bond] = {}
        self.MolsysInfo = MolsysInfoclass()
        self.CoeffsInfo = CoeffsInfoclass()

        # Molecular system basic type information
        self.MolsysInfo.angle_type = temp.angle_type
        self.MolsysInfo.angles_num = temp.angles_num
        self.MolsysInfo.atom_type = temp.atom_type
        self.MolsysInfo.atoms_num = temp.atoms_num
        self.MolsysInfo.bond_type = temp.bond_type
        self.MolsysInfo.bonds_num = temp.bonds_num
        self.MolsysInfo.dihedral_type = temp.dihedral_type
        self.MolsysInfo.dihedrals_num = temp.dihedrals_num
        self.MolsysInfo.xhi = temp.xhi
        self.MolsysInfo.xlo = temp.xlo
        self.MolsysInfo.yhi = temp.yhi
        self.MolsysInfo.ylo = temp.ylo
        self.MolsysInfo.zhi = temp.zhi
        self.MolsysInfo.zlo = temp.zlo

        # Atoms, Bonds, etc.

        self.Atoms: dict[int, Polymer.atom] = temp.Atoms
        self.Bonds: dict[int, Polymer.bond] = temp.Bonds
        self.Angles: dict[int, Polymer.angle] = temp.Angles
        self.Dihedrals: dict[int, Polymer.dihedral] = temp.Dihedrals

        # On first read, file ordering may be messy, need to reorganize

        # something else
        if hasattr(temp, "Pair_Coeffs"):
            self.CoeffsInfo.Pair_Coeffs = temp.Pair_Coeffs
        if hasattr(temp, "Angle_Coeffs"):
            self.CoeffsInfo.Angle_Coeffs = temp.Angle_Coeffs

        if hasattr(temp, "Bond_Coeffs"):
            self.CoeffsInfo.Bond_Coeffs = temp.Bond_Coeffs
        if hasattr(temp, "Dihedral_Coeffs"):
            self.CoeffsInfo.Dihedral_Coeffs = temp.Dihedral_Coeffs
        if hasattr(temp, "PairIJ_Coeffs"):
            self.CoeffsInfo.PairIJ_Coeffs = temp.PairIJ_Coeffs
        if hasattr(temp, "Masses"):
            self.CoeffsInfo.Masses = temp.Masses

        # Check Impropers separately
        if hasattr(temp, "Impropers"):
            self.MolsysInfo.improper_type = temp.improper_type
            self.MolsysInfo.impropers_num = temp.impropers_num
            self.CoeffsInfo.Improper_Coeffs = temp.Improper_Coeffs
            self.Impropers = temp.Impropers

        del temp

        # Read check, verify read count matches declared count in file
        if (
            (self.MolsysInfo.impropers_num != len(self.Impropers))
            or (self.MolsysInfo.dihedrals_num != self.Dihedrals.__len__())
            or (self.MolsysInfo.bonds_num != self.Bonds.__len__())
            or (self.MolsysInfo.angles_num != self.Angles.__len__())
            or (self.MolsysInfo.atoms_num != self.Atoms.__len__())
        ):
            raise ValueError("The number of read objects does not match the declared number in the file!")

        self.setMolList()

    def getAtom(self, atom_str_or_index: int|str) -> Polymer.atom:
        """Get Atom object

        Args:
            atom_str_or_index (int | str): Atom index or type

        Raises:
            KeyError: Atom index not found
            KeyError: Atom type not found
            ValueError: Parameter must be int or str

        Returns:
            Polymer.atom: Return corresponding Atom object
        """
        if type(atom_str_or_index) == int:
            if 1 <= atom_str_or_index <= self.MolsysInfo.atoms_num:
                return self.Atoms[atom_str_or_index]
            else:
                raise KeyError(f"No such Atom which index is {atom_str_or_index}")
        elif type(atom_str_or_index) == str:
            for atom in self.Atoms.values():
                if atom.abbreviation == atom_str_or_index:
                    return atom
            raise KeyError(f"No such Atom which type is {atom_str_or_index}")
        else:
            raise ValueError("arg must be int or str")

    def getBond(self, index) -> Polymer.bond:
        if 1 <= index <= self.MolsysInfo.bonds_num:
            return self.Bonds[index]
        else:
            raise KeyError(f"No such Bond which index is {index}")

    def getAngle(self, index) -> Polymer.angle:
        if 1 <= index <= self.MolsysInfo.angles_num:
            return self.Angles[index]
        else:
            raise KeyError(f"No such Angle which index is {index}")

    def getDihedral(self, index) -> Polymer.dihedral:
        if 1 <= index <= self.MolsysInfo.dihedrals_num:
            return self.Dihedrals[index]
        else:
            raise KeyError(f"No such Dihedral which index is {index}")

    def getImproper(self, index) -> Polymer.improper:
        if 1 <= index <= self.MolsysInfo.impropers_num:
            return self.Impropers[index]
        else:
            raise KeyError(f"No such Improper which index is {index}")
        
    def groupAtoms(self, group: str|int) -> list[Polymer.atom]:
        if type(group) == int:
            return [atom for atom in self.Atoms.values() if atom.type == group]
        elif type(group) == str:
            return [atom for atom in self.Atoms.values() if atom.abbreviation == group]
        else:
            raise ValueError("group must be str or int")

    def setMolList(self) -> None:
        """Set molecule list"""
        self.Mols: dict[int, list[Polymer.atom]] = {}
        for k, v in self.Atoms.items():
            mol_id: int = v.mol
            if mol_id not in self.Mols:
                self.Mols[mol_id] = []
            self.Mols[mol_id].append(self.Atoms[k])

    def update(self) -> None:
        self.MolsysInfo.angles_num = len(self.Angles)
        self.MolsysInfo.atoms_num = len(self.Atoms)
        self.MolsysInfo.bonds_num = len(self.Bonds)
        self.MolsysInfo.dihedrals_num = len(self.Dihedrals)
        self.MolsysInfo.impropers_num = len(self.Impropers)

    def writeNewLmpData(self, PATH) -> None:
        self.update()

        writedata: list[str] = [
            "Writed by IPS",
            " ",
        ]
        writedata += [
            f'{getattr(self.MolsysInfo, attr)} {attr.split("_")[0]}'
            for attr in [
                "atoms_num",
                "bonds_num",
                "angles_num",
                "dihedrals_num",
                "impropers_num",
            ]
        ]
        writedata += [""]
        writedata += [
            f'{getattr(self.MolsysInfo, attr)} {attr.split("_")[0]} types'
            for attr in [
                "atom_type",
                "bond_type",
                "angle_type",
                "dihedral_type",
                "improper_type",
            ]
        ]
        writedata += [""]
        writedata += [
            f"{getattr(self.MolsysInfo, attr)} {getattr(self.MolsysInfo, attr.replace('lo', 'hi'))} {attr.split('lo')[0]}lo {attr.split('lo')[0]}hi"
            for attr in ["xlo", "ylo", "zlo"]
        ]
        writedata += [""]

        for attr in [
            "Masses",
            "Pair_Coeffs",
            "PairIJ_Coeffs",
            "Bond_Coeffs",
            "Angle_Coeffs",
            "Dihedral_Coeffs",
            "Improper_Coeffs",
        ]:
            if getattr(self.CoeffsInfo, attr, None):
                writedata += [""]
                writedata += [attr.replace("_", " "), ""]
                writedata += [
                    " ".join(map(str, x)) for x in getattr(self.CoeffsInfo, attr)
                ]
                writedata += [""]

        for attr in ["Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"]:
            writedata += [attr, ""]
            if attr == "Atoms":
                writedata += [
                    f"{x.index} {x.mol} {x.type} {x.charge} {x.x} {x.y} {x.z}"
                    for k, x in getattr(self, attr).items()
                ]
            else:
                writedata += [
                    f"{x.index} {x.type} {' '.join(map(str, x.atoms_num))}"
                    for k, x in getattr(self, attr).items()
                ]
            writedata += [""]

        with open(PATH, "w", encoding="utf-8") as w:
            w.writelines(map(lambda s: s + "\n", writedata))
        # print(f"File writed! Write in {PATH}")


class readlammps(object):
    def __init__(self, lammpsfile) -> None:
        self.Atoms = {}
        self.Labels = []
        self.readLammpsfile(lammpsfile)

    def set_atoms(
        self,
        atomsList: list[Polymer.multi],
        bondslist: list[Polymer.multi],
        Angles: list[Polymer.multi],
        Dihedrals: list[Polymer.multi],
        Impropers: list[Polymer.multi] = None,
    ) -> None:
        self.set_multi_properties(bondslist, "bonds", "bonded")
        self.set_multi_properties(Angles, "angles")
        self.set_multi_properties(Dihedrals, "diheds")
        if Impropers:
            self.set_multi_properties(Impropers, "improps")

    def set_multi_properties(
        self, multi_list: dict[Polymer.multi], property_name, bonded_property=None
    ):
        for k, v in multi_list.items():
            for atom in v.atoms:
                getattr(atom, property_name).append(multi_list[k])
                if bonded_property:
                    getattr(atom, bonded_property).extend(
                        [i for i in v.atoms if i.index != atom.index]
                    )

    def set_multi(
        self,
        bondList: list,
        AnglesList: list,
        DihedralsList: list,
        ImpropersList: list = None,
    ):
        self.Bonds = {
            k: Polymer.bond(v[:2] + [self.Atoms[i] for i in v[2:]])
            for k, v in bondList.items()
        }
        self.Angles = {
            k: Polymer.angle(v[:2] + [self.Atoms[i] for i in v[2:]])
            for k, v in AnglesList.items()
        }
        self.Dihedrals = {
            k: Polymer.dihedral(v[:2] + [self.Atoms[i] for i in v[2:]])
            for k, v in DihedralsList.items()
        }
        if ImpropersList:
            self.Impropers = {
                k: Polymer.improper(v[:2] + [self.Atoms[i] for i in v[2:]])
                for k, v in ImpropersList.items()
            }

    def readLammpsfile(self, filename):
        with open(filename, "r") as file:
            data = file.readlines()
        for lineindex in range(len(data)):
            if len(data[lineindex]) != 0 and data[lineindex][-1] == "\n":
                data[lineindex] = basicfunc.Dealine(data[lineindex])

        data = [x for x in data if len(x) != 0]

        for lineindex, line in enumerate(data):
            lenline = len(line)

            if lenline > 4:
                continue
            elif lenline == 1:
                if line[0] in [
                    "Atoms",
                    "Bonds",
                    "Angles",
                    "Dihedrals",
                    "Impropers",
                    "Masses",
                    "Velocities",
                ]:
                    self.Labels.append([lineindex, line[0]])
            elif lenline == 2:
                if line[1] in [
                    "atoms",
                    "bonds",
                    "angles",
                    "dihedrals",
                    "bonds",
                    "impropers",
                ]:
                    setattr(self, f"{line[1]}_num", line[0])
                elif (
                    line[0]
                    in [
                        "Pair",
                        "Bond",
                        "Angle",
                        "BondBond",
                        "BondAngle",
                        "Dihedral",
                        "MiddleBondTorsion",
                        "EndBondTorsion",
                        "AngleTorsion",
                        "AngleAngleTorsion",
                        "BondBond13",
                        "Improper",
                        "AngleAngle",
                        "PairIJ",
                    ]
                    and line[1] == "Coeffs"
                ):
                    self.Labels.append([lineindex, f"{line[0]}_Coeffs"])
            elif lenline == 3:
                if line[1] in ["atom", "bond", "angle", "dihedral", "improper"]:
                    setattr(self, f"{line[1]}_type", int(line[0]))
            elif lenline == 4:
                if line[2] in ["xlo", "ylo", "zlo"]:
                    setattr(self, line[2], line[0])
                    setattr(self, line[2].replace("lo", "hi"), line[1])

        self.Labels.sort()
        if self.Labels:
            """
            Initial data setup for atoms, bonds, angles, dihedrals, and impropers
            """
            for x in range(len(self.Labels) - 1):
                if self.Labels[x][1] in [
                    "Atoms",
                    "Velocities",
                    "Bonds",
                    "Angles",
                    "Dihedrals",
                    "Impropers",
                ]:
                    setattr(
                        self,
                        self.Labels[x][1],
                        {
                            item[0]: item
                            for i, item in enumerate(
                                data[self.Labels[x][0] + 1 : self.Labels[x + 1][0]],
                                start=1,
                            )
                        },
                    )
                else:
                    setattr(
                        self,
                        self.Labels[x][1],
                        data[self.Labels[x][0] + 1 : self.Labels[x + 1][0]],
                    )
            if self.Labels[-1][1] in [
                "Atoms",
                "Velocities",
                "Bonds",
                "Angles",
                "Dihedrals",
                "Impropers",
            ]:
                setattr(
                    self,
                    self.Labels[-1][1],
                    {
                        item[0]: item
                        for i, item in enumerate(
                            data[self.Labels[-1][0] + 1 :], start=1
                        )
                    },
                )
            else:
                setattr(self, self.Labels[-1][1], data[self.Labels[-1][0] + 1 :])
        # Safety check, some files don't have Impropers

        properties = ["Bonds", "Angles", "Dihedrals", "Impropers"]
        existing_properties = [prop for prop in properties if hasattr(self, prop)]

        if len(existing_properties) < 3:
            raise basicfunc.DefinedError("Bonds or Angles not defined!")

        # Setup Atoms
        self.Atoms = {k: Polymer.atom(v) for k, v in self.Atoms.items()}

        # Setup Bonds, Angles, Dihedrals, Impropers
        self.set_multi(*[getattr(self, prop) for prop in existing_properties])

        # Setup bond, angle, dihedral, improper properties for Atoms
        self.set_atoms(
            self.Atoms, *[getattr(self, prop) for prop in existing_properties]
        )


