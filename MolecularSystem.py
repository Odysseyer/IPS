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

from readlammps import lammps
import Polymer
import readConfig
import numpy as np
from basicfunc import gift_wrapping, find_valid_placement


class MolecularSystem:
    """Molecular system class"""

    def __init__(self, lmpdatafile, config: dict) -> None:
        """Initialize molecular system

        Args:
            lmpdatafile (_type_): LAMMPS data file
            config (dict): Configuration file
        """

        lmpdata = lammps(lmpdatafile)

        self.types: readConfig.types = config["typedata"]
        self.instruction: readConfig.instruction = config["instruction"]
        self.lmp: lammps = lmpdata

        self.LoadAbbreviation()

        if self.instruction.equilibrium:
            self.upperPlate: Polymer.atom = self.lmp.getAtom(self.instruction.upperPlate)
            self.lowerPlate: Polymer.atom = self.lmp.getAtom(self.instruction.lowerPlate)

            self.upperSolvent = set(self.instruction.upperSolvent)
            self.upperMonomer: Polymer.Molecular = Polymer.Molecular()
            self.upperMonomer.readMononerfile(self.instruction.upperMonomer)

            self.lowerSolvent = set(self.instruction.lowerSolvent)
            self.lowerMonomer: Polymer.Molecular = Polymer.Molecular()
            self.lowerMonomer.readMononerfile(self.instruction.lowerMonomer)

    def LoadAbbreviation(self) -> None:
        """
        Load abbreviations in the system
        """
        for atom in self.lmp.Atoms.values():
            atom.findtypestr(self.types.atom_Type)

    def queryType(self, items: Polymer.multi) -> str:
        """Query types of bonds, angles, and dihedrals
        Args:
            items (Polymer.multi): Bond, angle, or dihedral object
        Raises:
            KeyError: Type not found
            KeyError: Type not found
            KeyError: Type not found
            TypeError: Wrong type
        Returns:
            _type_: Type
        """
        Reverse_Order: list[str] = items.abbreviation.split(",")
        Reverse_Order.reverse()
        Reverse_Order = ",".join(Reverse_Order)
        Forward_Order: str = items.abbreviation

        if isinstance(items, Polymer.bond):
            if Forward_Order in self.types.bond_Type:
                return self.types.bond_Type[Forward_Order]

            elif Reverse_Order in self.types.bond_Type:
                return self.types.bond_Type[Reverse_Order]

            else:
                raise KeyError(f"Not find this type {Forward_Order}")

        elif isinstance(items, Polymer.angle):

            if Forward_Order in self.types.angle_Type:
                return self.types.angle_Type[Forward_Order]

            elif Reverse_Order in self.types.angle_Type:
                return self.types.angle_Type[Reverse_Order]

            else:
                raise KeyError(f"Not find this type {Forward_Order}")
        elif isinstance(items, Polymer.dihedral):

            if Forward_Order in self.types.dihedral_Type:
                return self.types.dihedral_Type[Forward_Order]

            elif Reverse_Order in self.types.dihedral_Type:
                return self.types.dihedral_Type[Reverse_Order]

            else:
                raise KeyError(f"Not find this type {Forward_Order}")

        elif isinstance(items, Polymer.improper):
            if Forward_Order in self.types.improper_Type:
                return self.types.improper_Type[Forward_Order]
            else:
                raise KeyError(f"Not find this type {Forward_Order}")
        else:
            raise TypeError(
                f"Give Wrong item!, it is not a improper/dihedral/angle/bond"
            )

    def queryType_new(self, items: list[Polymer.atom], type) -> str:
        """Query types of bonds, angles, and dihedrals

        Args:
            items (list[Polymer.atom]): List of atoms
            type (_type_): Type

        Raises:
            KeyError: Type not found
            KeyError: Type not found
            KeyError: Type not found
            TypeError: Wrong type

        Returns:
            str: Type
        """

        tempOrder: list[str] = [atom.abbreviation for atom in items]
        Forward_Order: str = ",".join(tempOrder)
        Reverse_Order: list[str] = tempOrder
        Reverse_Order.reverse()
        Reverse_Order = ",".join(Reverse_Order)

        if type == "bond":
            if Forward_Order in self.types.bond_Type:
                return self.types.bond_Type[Forward_Order]
            elif Reverse_Order in self.types.bond_Type:
                return self.types.bond_Type[Reverse_Order]
            else:
                raise KeyError(f"Not find this type {Forward_Order}")

        elif type == "angle":
            if Forward_Order in self.types.angle_Type:
                return self.types.angle_Type[Forward_Order]
            elif Reverse_Order in self.types.angle_Type:
                return self.types.angle_Type[Reverse_Order]
            else:
                raise KeyError(f"Not find this type {Forward_Order}")

        elif type == "dihedral":
            if Forward_Order in self.types.dihedral_Type:
                return self.types.dihedral_Type[Forward_Order]
            elif Reverse_Order in self.types.dihedral_Type:
                return self.types.dihedral_Type[Reverse_Order]
            else:
                raise KeyError(f"Not find this type {Forward_Order}")
        else:
            raise TypeError(f"Give Wrong item!, it is not a dihedral/angle/bond")

    def queryType_Improper(self, improper: list[Polymer.atom]) -> str:
        """Query improper angle types

        Args:
            improper (list[Polymer.atom]): Atoms of improper angle

        Raises:
            KeyError: If improper angle type is not found

        Returns:
            _type_: Improper angle type
        """
        c, a, b,  d = improper
        # Require central atom in first position, permute remaining atoms with six results
        permutation_order: list[list[Polymer.atom]] = [
            [c, a, b,  d],
            [c, b, a,  d],
            [c, d, b,  a],
            [c, b, d,  a],
            [c, a, d,  b],
            [c, d, a,  b],
        ]

        for order in permutation_order:
            ordstr: list[str] = [atom.abbreviation for atom in order]
            ordstr = ",".join(ordstr)
            if ordstr in self.types.improper_Type:
                return self.types.improper_Type[ordstr]
        raise KeyError(
                f"Not find this improper's type, {c.abbreviation},{a.abbreviation},{b.abbreviation},{d.abbreviation}"
            )

    def merge(self, other: Polymer.Molecular, direction: str) -> None:
        """
        Merge new molecule into the system

        Args:
            other (Polymer.Molecular): Monomer to merge
            direction (str): Merge direction, 'up' means add monomer from above, 'be' means add monomer from below

        """

        if type(other) != Polymer.Molecular:
            raise TypeError("The other object must be Monomer object")

        def findVacancies(plateZ: int) -> np.array:
            """
            Find vacancies in current system
            Args:
                `plateZ` (int): Z coordinate of metal plate
            Returns:
                np.ndarray: Coordinates of vacancies
            """
            filtertypes=list(self.upperSolvent)+list(self.lowerSolvent)+list([self.types.atom_Type["be"],self.types.atom_Type["up"]])
            circles = np.array(
                [
                    (atom.x, atom.y)
                    for atom in self.lmp.Atoms.values()
                    if (plateZ - 38 <= atom.z <= plateZ + 38) and (atom.type not in filtertypes)
                ]
            )
            # Calculate outermost points
            points = np.array([(atom.x, atom.y) for atom in other.atoms.values()])
            hull = gift_wrapping(points)

            # Try to find valid placement position
            valid_points = find_valid_placement(
                hull, points, circles, 0.88, self.lmp.MolsysInfo.xhi, max_attempts=2000
            )

            return valid_points

        if direction == "up":
            plateZ = self.lmp.getAtom("up")
            valid_points = findVacancies(plateZ.z)
            for index, atom in other.atoms.items():
                atom.x, atom.y, atom.z = (
                    valid_points[index-1][0].item(),
                    valid_points[index-1][1].item(),
                    plateZ.z - 1,
                )  # Convert to Python native float

        elif direction == "be":
            plateZ = self.lmp.getAtom("be")
            valid_points = findVacancies(plateZ.z)
            for index, atom in other.atoms.items():
                atom.x, atom.y, atom.z = (
                    valid_points[index-1][0].item(),
                    valid_points[index-1][1].item(),
                    plateZ.z + 1,
                )

        else:
            raise ValueError("The direction must be up or down")

        # Update system information
        other.update(
            self.lmp.MolsysInfo.atoms_num,
            self.lmp.MolsysInfo.bonds_num,
            self.lmp.MolsysInfo.angles_num,
            self.lmp.MolsysInfo.dihedrals_num,
            self.lmp.MolsysInfo.impropers_num,
            self.lmp.Mols.__len__()
        )

        self.lmp.Atoms.update(other.atoms)
        self.lmp.Bonds.update(other.bonds)
        self.lmp.Angles.update(other.angles)
        self.lmp.Dihedrals.update(other.dihedrals)
        self.lmp.Impropers.update(other.impropers)
        self.lmp.Mols[self.lmp.Mols.__len__()+1]=[a for a in other.atoms.values()]

        self.lmp.MolsysInfo.atoms_num = len(self.lmp.Atoms)
        self.lmp.MolsysInfo.bonds_num = len(self.lmp.Bonds)
        self.lmp.MolsysInfo.angles_num = len(self.lmp.Angles)
        self.lmp.MolsysInfo.dihedrals_num = len(self.lmp.Dihedrals)
        self.lmp.MolsysInfo.impropers_num = len(self.lmp.Impropers)

    def deleteAtom(self, atom: Polymer.atom) -> None:
        """
        Delete all related objects of deleted atoms from the atomic system
        """
        
        for bonded_atom in atom.bonded:
            bonded_atom.bonded.remove(atom)

        for bond in atom.bonds:
            for atom_in_bond in bond.atoms:
                if atom_in_bond != atom:
                    atom_in_bond.bonds.remove(bond)
            self.lmp.Bonds.pop(bond.index)
            

        for angle in atom.angles:
            for atom_in_angle in angle.atoms:
                if atom_in_angle != atom:
                    atom_in_angle.angles.remove(angle)
            self.lmp.Angles.pop(angle.index)

        for dihed in atom.diheds:
            for atom_in_dihed in dihed.atoms:
                if atom_in_dihed != atom:
                    atom_in_dihed.diheds.remove(dihed)
            self.lmp.Dihedrals.pop(dihed.index)

        for improp in atom.improps:
            for atom_in_improp in improp.atoms:
                if atom_in_improp != atom:
                    atom_in_improp.improps.remove(improp)
            self.lmp.Impropers.pop(improp.index)
            
        self.lmp.Atoms.pop(atom.index)

    def equilibrium(self) -> bool:
        """Analyze system equilibrium

        Args:
            `upperMonomer` (Polymer.Molecular): Upper monomer
            `upperSolvent` (list[int]): Upper solvent
            `upperFactor` (int): Upper monomer-solvent ratio factor
            `lowerMonomer` (Polymer.Molecular): Lower monomer
            `lowerSolvent` (list[int]): Lower solvent
            `lowerFactor` (int): Lower monomer-solvent ratio factor
        """
        NEEDEQUILIBRIUM = False
        if self.analysis(
            self.upperMonomer.type, self.upperSolvent, self.instruction.upperFactor
        ):  # If monomer needs to be added on upper side
            self.merge(self.upperMonomer, "up")  # Merge monomer to upper side
            # Optimizable, pass copy of monomer instead of monomer itself
            self.movePlate("up", 0.065)  # Move upper metal plate, 0.065 is empirical value
            NEEDEQUILIBRIUM  = True
            print(f'[INFO]\tAdd A Monomer to the upper side.')

        if self.analysis(
            self.lowerMonomer.type, self.lowerSolvent, self.instruction.lowerFactor
        ):  # If monomer needs to be added on lower side
            self.merge(self.lowerMonomer, "be")  # Merge monomer to lower side
            self.movePlate("be", 0.065)  # Move lower metal plate, 0.065 is empirical value
            NEEDEQUILIBRIUM = True
            print(f'[INFO]\tAdd A Monomer to the lower side.')
            
        return NEEDEQUILIBRIUM

    def analysis(
        self, MonomerType: set[int], SolventType: set[int], factor: int
    ) -> bool:
        """Analyze whether monomers need to be added to the system

        Args:
            `MonomerType` (set[int]): Monomer type
            `SolventType` (set[int]): Solvent type
            `factor` (int): Monomer-solvent ratio factor

        Returns:
            `bool`: Whether monomers need to be added
        """
        def tran2strtype(typenum):
            strtype=[self.types.atom_Type[num] for num in typenum]
            return "-".join(strtype)
        
        def print_info(monomer_counts, solvent_counts, factor, equilibrium):
            product = monomer_counts * factor
            status = "Equilibrium" if equilibrium else "Not Equilibrium"
            print(f'[INFO]\tMonomer [{tran2strtype(MonomerType)}]: {monomer_counts}\n'
                f'\tSolvent [{tran2strtype(SolventType)}]: {solvent_counts}\n'
                f'\tFactor: {factor}\n'
                f'\tProduct (Monomer Counts * Factor): {product}\n'
                f'\tStatus: {status}\n')

        if self.instruction.link.__len__()==2:
            linkatoms: list[float]=[atom.z for atom in self.lmp.groupAtoms(self.instruction.link[0][0][1])+self.lmp.groupAtoms(self.instruction.link[0][1][1])+self.lmp.groupAtoms(self.instruction.link[1][0][1])+self.lmp.groupAtoms(self.instruction.link[1][1][1])]
        else:
            linkatoms: list[float]=[atom.z for atom in self.lmp.groupAtoms(self.instruction.link[0][0][1])+self.lmp.groupAtoms(self.instruction.link[0][1][1])]
        
        def find_main_peak_boundaries(linkatoms, bin_width=2, threshold_ratio=0.1, EXPANDINGCONSTANTS=10):
            # Create histogram
            hist, bin_edges = np.histogram(linkatoms, bins=np.arange(min(linkatoms), max(linkatoms) + bin_width, bin_width))

            # Find highest peak
            peak_index = np.argmax(hist)
            peak_value = bin_edges[peak_index]

            # Set threshold
            threshold = max(hist) * threshold_ratio

            # Search left for boundary
            left_boundary = peak_value
            for i in range(peak_index, 0, -1):
                if hist[i] < threshold:
                    left_boundary = bin_edges[i]
                    break

            # Search right for boundary
            right_boundary = peak_value
            for i in range(peak_index, len(hist)):
                if hist[i] < threshold:
                    right_boundary = bin_edges[i]
                    break
            
            return left_boundary-EXPANDINGCONSTANTS, right_boundary+EXPANDINGCONSTANTS
        
        left, right = find_main_peak_boundaries(linkatoms)
        
        SearchSpace: list[int]=[atom.mol for atom in self.lmp.Atoms.values() if not (left <= atom.z <= right)]
        SearchSpace = set(SearchSpace)
        
        MonomerCounts = 0
        SolventCounts = 0
        
        for mol in SearchSpace:
            iterType = set([a.type for a in self.lmp.Mols[mol]])
            if iterType == MonomerType:
                MonomerCounts += 1
            elif iterType == set(SolventType):
                SolventCounts += 1
        
        if MonomerCounts * factor < SolventCounts:
            print_info(MonomerCounts,SolventCounts,factor,False)
            return True
        else:
            print_info(MonomerCounts,SolventCounts,factor,True)
            return False

    def movePlate(self, direction: str, step: int) -> None:
        """
        Move metal plate

        Args:
            `plate` (str): Metal plate name
            `direction` (str): Move direction
            `step` (int): Move step size
        """
        if direction == "up":
            for atom in self.lmp.Atoms.values():
                if atom.abbreviation == self.instruction.upperPlate:
                    atom.z += step
            self.lmp.MolsysInfo.zhi += step
        elif direction == "be":
            for atom in self.lmp.Atoms.values():
                if atom.abbreviation == self.instruction.lowerPlate:
                    atom.z -= step
            self.lmp.MolsysInfo.zlo -= step
        else:
            raise ValueError("The direction must be `up` or `be`")


