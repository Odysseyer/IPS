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

import itertools
from MolecularSystem import MolecularSystem
from Polymer import atom, bond, angle, dihedral, improper
import numpy as np # type: ignore



class InterfaceReaction(object):

    def __init__(self, lmpdatafile, config: dict, outfile=None) -> None:
        self.MolecularSystem = MolecularSystem(
            lmpdatafile, config)
        if outfile:
            self.outfile = outfile


    def finalDeal(self):

        # flag1 = False
        # flag2 = False
        # for atom in self.MolecularSystem.lmp.Atoms.values():
        #     if atom.abbreviation == self.MolecularSystem.instruction.link[0][0][0]:
        #         atom.charge -= self.MolecularSystem.instruction.charge[0]
        #         flag1 = True
        #     elif atom.abbreviation == self.MolecularSystem.instruction.link[0][1][0]:
        #         atom.charge -= self.MolecularSystem.instruction.charge[1]
        #         flag2 = True

        # if (not flag2) or (not flag1):
        #     raise KeyError(
        #         f"Do not find link atom which type is {self.MolecularSystem.instruction.link[0][0]}")
        

        self.MolecularSystem.lmp.writeNewLmpData(self.outfile)

    def initCharges(self) -> None:

        # flag1 = False
        # flag2 = False
        # # Assign given positive/negative charges to first atoms of link types 1 and 2
        # for atom in self.MolecularSystem.lmp.Atoms.values():
        #     if atom.abbreviation == self.MolecularSystem.instruction.link[0][0][0]:
        #         atom.charge += self.MolecularSystem.instruction.charge[0]
        #         flag1 = True
        #     elif atom.abbreviation == self.MolecularSystem.instruction.link[0][1][0]:
        #         atom.charge += self.MolecularSystem.instruction.charge[1]
        #         flag2 = True

        # if (not flag2) or (not flag1):
        #     raise KeyError(
        #         f"Do not find link atom which type is {self.MolecularSystem.instruction.link[0]}")
        self.MolecularSystem.lmp.writeNewLmpData(self.outfile)
        # print("[INFO]\tInitialization complete, Charge loaded!")
        return 0

    def findPair(self,link) -> list[atom]:
        """
        Find atom pairs that meet bonding criteria
        """

        # linkAtoms_A/B are sets of atoms matching link atom type numbers from .in file
        # Get connecting atom numeric types, find connecting atoms in system
        linkAtoms_A: list[atom] = self.MolecularSystem.lmp.groupAtoms(link[0][0])
        linkAtoms_B: list[atom] = self.MolecularSystem.lmp.groupAtoms(link[1][0])

        closest = 0

        # Check bonding criteria for all paired connecting atoms
        # x,y are atoms in linkAtoms_A/B respectively
        try:
            
            for x, y in itertools.product(linkAtoms_A, linkAtoms_B):

                # Intra check - check internal connecting atoms
                if self.MolecularSystem.instruction.intra:
                    # Check if given pair is within N bonds
                    if (self.checkbond(x, y)):
                        print('uuuu')
                        continue
                elif (x.mol == y.mol):
                    continue

                # Cutoff check - distance cutoff check
                sep: float = self.getSep(x, y)

                if sep > self.MolecularSystem.instruction.cutoff:
                    continue
                # Alignment check
                if self.MolecularSystem.instruction.align:
                    # Perform alignment checks for given pair
                    if (not self.aligned()):
                        continue
                # If this is the closest pair of connecting atoms, save this pair

                if (closest == 0 or sep < closest):
                    closest: float = sep
                    pair: list[atom] = [x, y]
        except Exception as e:
            print(e)

        # Return closest pair of connecting atoms
        if closest == 0:
            return False, 999
        else:
            print(
                f"\033[36m[INFO]\tPair:{pair} {closest:.4f} Å\033[0m")
            return pair,closest

    def checkbond(self, atom1: atom, atom2: atom) -> bool:
        """
        Check if atom2 exists in atom1's bond structure
        """
        # ----Breadth-first search on N-ary tree-----
        next_node = []
        node: list[atom] = [atom1]
        counts = 0
        while node:
            if counts >= 6:
                return False
            for x in node:
                if atom2 in x.bonded:
                    return True
                for y in x.bonded:
                    next_node.append(y)
            node = next_node
            next_node = []
            counts += 1

    def getSep(self, atom1: atom, atom2: atom) -> float:
        """
        Calculate distance between two given atoms
        """

        # Use vector operations here
        # Calculate difference between two atom positions
        position = atom1.position-atom2.position

        # length is the box size
        length = np.array([self.MolecularSystem.lmp.MolsysInfo.xhi-self.MolecularSystem.lmp.MolsysInfo.xlo,
                          self.MolecularSystem.lmp.MolsysInfo.yhi-self.MolecularSystem.lmp.MolsysInfo.ylo, self.MolecularSystem.lmp.MolsysInfo.zhi-self.MolecularSystem.lmp.MolsysInfo.zlo])

        # Move coordinates to center
        position = position - length*np.floor(position/length+0.5)

        # Return vector magnitude, which is the distance
        return np.linalg.norm(position)

    def get_connect(self, pair: list[atom], graph) -> list[None]:
        """
        Get atoms described around the bonding atoms in the command
        """
        queue: list[int] = [0, 1]
        result_atoms: list[None] = [None] * (len(graph))

        result_atoms[0] = pair[0]
        result_atoms[1] = pair[1]

        while queue:
            n: int = queue.pop(0)
            connect = graph[n]['bonded']
            bonded = result_atoms[n].bonded

            for atom in bonded: 
                for atomNum in connect:
                    if atom.type == graph[atomNum]['type']:
                        if result_atoms[atomNum] is not None:
                            raise ValueError(
                                "Atom connectivity in input script is not unique.")
                        result_atoms[atomNum] = atom
                        if 'bonded' in graph[atomNum]:
                            queue.append(atomNum)

        return result_atoms

    def aligned(self, pair) -> None | bool:
        """
        Check if given bond meets angle requirements
        """
        atoms = self.get_connect(pair)

        for item in self.MolecularSystem.instruction.align:
            atomA = [atoms[index].postion for index in item[0]]
            atomB = [atoms[index] for index in item[1]]

            result_angle = self.angle_between_vectors_or_planes(atomA, atomB)

            # This code has a security vulnerability
            if eval(item[-1].replace("angle", result_angle)):
                return True
            else:
                return False

    def getExtraBonds(self, pair):
        """
        Get extra bonds
        """
        atoms = self.get_connect(pair)

        extraBonds = []

        for bond in range(self.MolecularSystem.instruction.bond):
            atomA = atoms[bond[0]]
            atomB = atoms[bond[1]]
            extraBonds.append(atomA, atomB)
        return extraBonds

    def angle_between_vectors_or_planes(a: list, b: list) -> np.float64:
        """
        Calculate angle between given vectors or planes, including vector-vector, vector-plane, and plane-plane angles.
        Parameters a and b are vectors or planes for angle calculation, formatted as lists with coordinate points.
        Returns the angle value in degrees as float.
        """
        def fit_plane(vector: np.array) -> np.array:
            """
            Use PCA on given coordinate points: after representing coordinate points as a matrix,
            obtain the fitted plane normal vector by solving for the eigenvectors of the covariance matrix,
            and return the plane normal vector.
            Parameter vector is the coordinate point set for plane fitting, type is 2D numpy array, each row is a 3D coordinate point.
            Returns the fitted plane normal vector as 1D numpy array.
            """
            # Calculate covariance matrix and eigenvectors
            covariance_matrix = np.cov(vector.T)
            eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)

            # Select eigenvector corresponding to smallest eigenvalue as plane normal
            plane_normal = eigenvectors[:, np.argmin(eigenvalues)]

            return plane_normal

        if len(a) == 2:
            vectorA = a[0]-a[1]
        else:
            vectorA = fit_plane(np.vstack(a))

        if len(b) == 2:
            vectorB = b[0]-b[1]
        else:
            vectorB = fit_plane(np.vstack(b))

        # Calculate angle between two normal vectors
        cos_theta = np.dot(vectorA, vectorB) / \
            (np.linalg.norm(vectorA) * np.linalg.norm(vectorB))
        angle = np.arccos(cos_theta) * 180 / np.pi

        return angle

    def deleterealtions(self, pair: list[atom]) -> None:
        """
        Delete specific atoms connected to bonding atoms
        """
        for delAtom in pair[0].bonded:
            if delAtom.abbreviation in self.MolecularSystem.instruction.deletelink:
                print(f"[INFO]\tDelete Atom {delAtom.abbreviation} {delAtom.index}")
                self.MolecularSystem.deleteAtom(delAtom)
                break
            
        for delAtom in pair[1].bonded:
            if delAtom.abbreviation in self.MolecularSystem.instruction.deletelink:
                print(f"[INFO]\tDelete Atom {delAtom.abbreviation} {delAtom.index}")
                self.MolecularSystem.deleteAtom(delAtom)
                break

    def Update(self, pair: list[atom],link) -> None:
        """
        Update related object information after bonding to simulate bond formation
        """

        NewBonds: list = [pair]
        NewAngles: list = []
        NewDiheds: list = []
        NewImprops: list = []
        
        self.deleterealtions(pair)

        # Read atom types dictionary from types.txt, convert link commands from in file to corresponding numbers

        AtomAType:int = self.MolecularSystem.types.atom_Type[link[0][0]]
        newAtomAType:int = self.MolecularSystem.types.atom_Type[link[0][1]]
        AtomBType:int = self.MolecularSystem.types.atom_Type[link[1][0]]
        newAtomBType:int = self.MolecularSystem.types.atom_Type[link[1][1]]


        atomA:atom = pair[0]
        atomB:atom = pair[1]
        atomA.bonded.append(pair[1])
        atomB.bonded.append(pair[0])

        # Determine if charge modification command exists and atom types meet charge addition/subtraction conditions
        if self.MolecularSystem.instruction.charge and atomA.type == AtomAType and atomB.type == AtomBType:
            atomA.charge -= self.MolecularSystem.instruction.charge[0]
            atomB.charge -= self.MolecularSystem.instruction.charge[1]

        # Update atom types
        if atomA.type == AtomAType:
            atomA.type = newAtomAType

        if atomB.type == AtomBType:
            atomB.type = newAtomBType

        # Update changed atom's typestr
        atomA.findtypestr(self.MolecularSystem.types.atom_Type)
        atomB.findtypestr(self.MolecularSystem.types.atom_Type)

        # Update affected bonded terms
        # Bond
        update:list[atom] = list(set(atomA.bonds+atomB.bonds))  # Remove duplicate elements
        for y in update:
            # Find bond in update list 
            self.MolecularSystem.lmp.getBond(y.index).type = self.MolecularSystem.queryType(y) # Update bond type
            
        # Angles
        update = list(set(atomA.angles+atomB.angles))
        for y in update:
            self.MolecularSystem.lmp.getAngle(y.index).type = self.MolecularSystem.queryType(y) # Update angle type
        # Dihedrals
        update = list(set(atomA.diheds+atomB.diheds))
        for y in update:
            self.MolecularSystem.lmp.getDihedral(y.index).type = self.MolecularSystem.queryType(y) # Update dihedral type
        # Impropers
        update = list(set(atomA.improps+atomB.improps))
        for y in update:
            self.MolecularSystem.lmp.getImproper(y.index).type = self.MolecularSystem.queryType(y) # Update improper type

        # Add new bonded terms
        # Generate new bond, angle, dihedral, improper
        """
        
        """
        if self.MolecularSystem.types.angle_Type == 0 and self.MolecularSystem.types.dihedral_Type == 0 and self.MolecularSystem.types.improper_Type == 0:
            raise ValueError("No angle, dihedral, or improper types found")

        bondedA: list[atom] = atomA.bonded
        bondedB: list[atom] = atomB.bonded
        
        # Add new angle, dihedral, composed of pair atoms and atoms connected to reaction atoms
        for atomC in atomB.bonded: # First search in B
            if (atomC == atomA):
                continue
            NewAngles.append([atomA, atomB, atomC])
            for atomD in atomC.bonded: # Search for atoms in C to form dihedral
                if atomD == atomB:
                    continue
                NewDiheds.append([atomA, atomB, atomC, atomD]) # A-B-X-X
                
        for atomC in bondedA: # Search in A
            if (atomC == atomB):
                continue
            NewAngles.append([atomB, atomA, atomC])

            for atomD in atomC.bonded: # Search for atoms in C to form dihedral
                if atomD == atomA:
                    continue
                NewDiheds.append([atomB, atomA, atomC, atomD]) # B-A-X-X

            for atomD in atomB.bonded: # Two reaction atoms in middle, only search once since X-A-B-X == X-B-A-X
                if atomD == atomA:
                    continue
                NewDiheds.append([atomC, atomA, atomB, atomD]) # X-A-B-X == X-B-A-X
        """
        Define improper torsion angles, here the first position is the central atom,
        others are bonded atoms of the central atom, only generated when central atom is bonded to three atoms
        """
        
        if self.MolecularSystem.types.improper_Type and len(atomA.bonded) == 3:
            NewImprops.append([atomA,atomA.bonded[0],atomA.bonded[1],  atomA.bonded[2]])
        if self.MolecularSystem.types.improper_Type and len(atomB.bonded) == 3:
            NewImprops.append([atomB,atomB.bonded[0],atomB.bonded[1],  atomB.bonded[2]])

        """Add new bonded terms to the system
        """

        for x in NewBonds:
            self.MolecularSystem.lmp.MolsysInfo.bonds_num += 1
            newType = self.MolecularSystem.queryType_new(x, 'bond')
            temp = [self.MolecularSystem.lmp.MolsysInfo.bonds_num, newType]
            temp = temp+[y for y in x]
            self.MolecularSystem.lmp.Bonds[self.MolecularSystem.lmp.MolsysInfo.bonds_num]=bond(temp)

        for x in NewAngles:
            self.MolecularSystem.lmp.MolsysInfo.angles_num += 1
            newType = self.MolecularSystem.queryType_new(x, 'angle')
            temp = [self.MolecularSystem.lmp.MolsysInfo.angles_num, newType]
            temp = temp+[y for y in x]
            self.MolecularSystem.lmp.Angles[self.MolecularSystem.lmp.MolsysInfo.angles_num]=angle(temp)

        for x in NewDiheds:
            self.MolecularSystem.lmp.MolsysInfo.dihedrals_num += 1
            newType: str = self.MolecularSystem.queryType_new(x, 'dihedral')
            temp = [self.MolecularSystem.lmp.MolsysInfo.dihedrals_num, newType]
            temp = temp+[y for y in x]
            self.MolecularSystem.lmp.Dihedrals[self.MolecularSystem.lmp.MolsysInfo.dihedrals_num]=dihedral(temp)

        for x in NewImprops:
            self.MolecularSystem.lmp.MolsysInfo.impropers_num += 1
            newType = self.MolecularSystem.queryType_Improper(x)
            temp = [self.MolecularSystem.lmp.MolsysInfo.impropers_num, newType]
            temp = temp+[y for y in x]
            self.MolecularSystem.lmp.Impropers[self.MolecularSystem.lmp.MolsysInfo.impropers_num]=improper(temp)

        self.MolecularUpdate(pair)
        

    def MolecularUpdate(self, pair: list[atom]) -> None:
        """
        First, check if two atoms belong to the same molecule, if so, do nothing and exit.
        If two atoms belong to different molecules, connect the two molecules into one large molecule.
        First, calculate the current number of molecules to move the last molecule to the deleted position.
        Next, set the molecule with smaller number as target molecule, larger number as source molecule.
        Then, move all atoms from source molecule to target molecule.
        Then, move atoms from the last molecule to the source molecule.
        Finally, delete the last molecule and update molecule counter.
        """
        if pair[0].mol != pair[1].mol:
            if pair[0].mol > pair[1].mol:
                maxMol, minMol = pair[0].mol, pair[1].mol
            else:
                minMol, maxMol = pair[0].mol, pair[1].mol

            # Move all atoms from larger molecule to smaller molecule
            for atom in self.MolecularSystem.lmp.Mols[maxMol]:
                atom.mol = minMol
                self.MolecularSystem.lmp.Mols[minMol].append(atom)

            # Clear larger molecule and fill with the last molecule
            self.MolecularSystem.lmp.Mols[maxMol] = []
            self.MolecularSystem.lmp.Mols[maxMol] = self.MolecularSystem.lmp.Mols[maxMol] + \
                self.MolecularSystem.lmp.Mols.pop(
                    max(self.MolecularSystem.lmp.Mols.keys()))
            # Finally, modify added atoms' Mol to be consistent
            for atom in self.MolecularSystem.lmp.Mols[maxMol]:
                atom.mol = maxMol
