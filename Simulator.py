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

import shutil
import yaml
import readConfig
import os
import subprocess
import re
from InterfaceReaction import InterfaceReaction
from datetime import datetime

class membranepolymer():

    def __init__(self):
        self.head()
        self.readconfig()

        if not os.path.isdir(self.outPath):
            os.makedirs(self.outPath)
        if not os.path.isdir(self.lmprundir):
            os.makedirs(self.lmprundir)
        if not os.path.isdir("trajectory"):
            os.makedirs("trajectory")
        # Initialize charges
        if self.bondsNum == 0:
            self.Prepare()
            self.run_lammps(self.minfile)
            shutil.move(f'{self.lmprundir}/min.lmps',
                        f'{self.outPath}/bond_0_0.lmps')

        while (self.bondsNum < self.bondLimits): # 
            result: bool = self.form_bonds()
            if not result or self.bondsNum == self.bondLimits:
                break
            self.get_newestfile(self.outPath, self.lmprundir)
            self.run_lammps(self.minfile)
            shutil.move(f'{self.lmprundir}/min.lmps',
                        f'{self.outPath}/bond_{self.bondsNum}_0.lmps')
            shutil.move(f'{self.lmprundir}/log.lammps',f'trajectory/bond_{self.bondsNum}_log.lammps')

            
        # Final deal, remove charges
        self.get_newestfile(self.outPath, self.lmprundir)
        self.run_lammps(self.minfile)
        shutil.move(f'{self.lmprundir}/min.lmps',
                    f'bond_final.lmps')
        shutil.move(f'{self.lmprundir}/log.lammps',f'trajectory/bond_final_log.lammps')

        # self.remove_charges()

    def head(self) -> None:

        print('--------------------------------------------------------------------------------------------------------------------------')
        print('            ██╗███╗   ██╗████████╗███████╗██████╗ ███████╗ █████╗  ██████╗██╗ █████╗ ██╗                                    ')
        print('            ██║████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔════╝██╔══██╗██╔════╝██║██╔══██╗██║                                    ')
        print('            ██║██╔██╗ ██║   ██║   █████╗  ██████╔╝█████╗  ███████║██║     ██║███████║██║                                    ')
        print('            ██║██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██╔══╝  ██╔══██║██║     ██║██╔══██║██║                                    ')
        print('            ██║██║ ╚████║   ██║   ███████╗██║  ██║██║     ██║  ██║╚██████╗██║██║  ██║███████╗                               ')
        print('            ╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═╝╚══════╝                               ')
        print('                                                                                                                ')
        print('██████╗  ██████╗ ██╗  ██╗   ██╗███╗   ███╗███████╗██████╗ ██╗███████╗ █████╗ ████████╗██╗ ██████╗ ███╗   ██╗    ')
        print('██╔══██╗██╔═══██╗██║  ╚██╗ ██╔╝████╗ ████║██╔════╝██╔══██╗██║╚══███╔╝██╔══██╗╚══██╔══╝██║██╔═══██╗████╗  ██║    ')
        print('██████╔╝██║   ██║██║   ╚████╔╝ ██╔████╔██║█████╗  ██████╔╝██║  ███╔╝ ███████║   ██║   ██║██║   ██║██╔██╗ ██║    ')
        print('██╔═══╝ ██║   ██║██║    ╚██╔╝  ██║╚██╔╝██║██╔══╝  ██╔══██╗██║ ███╔╝  ██╔══██║   ██║   ██║██║   ██║██║╚██╗██║    ')
        print('██║     ╚██████╔╝███████╗██║   ██║ ╚═╝ ██║███████╗██║  ██║██║███████╗██║  ██║   ██║   ██║╚██████╔╝██║ ╚████║    ')
        print('╚═╝      ╚═════╝ ╚══════╝╚═╝   ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝╚═╝╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝    ')
        print('                                                                                                                ')
        print('                ███████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗ ██████╗ ██████╗                                        ')
        print('                ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██╔═══██╗██╔══██╗                                       ')
        print('                ███████╗██║██╔████╔██║██║   ██║██║     ███████║   ██║   ██║   ██║██████╔╝                                       ')
        print('                ╚════██║██║██║╚██╔╝██║██║   ██║██║     ██╔══██║   ██║   ██║   ██║██╔══██╗                                       ')
        print('                ███████║██║██║ ╚═╝ ██║╚██████╔╝███████╗██║  ██║   ██║   ╚██████╔╝██║  ██║                                       ')
        print('                ╚══════╝╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝                                       ')
        print('--------------------------------------------------------------------------------------------------------------------------')

    def form_bonds(self) -> bool:
        while True:
            print(f'-'*10)
            print(f'[INFO]\tBonds:{self.bondsNum}\tTryCounts:{self.tryNum}')
            Mem = InterfaceReaction(
                self.find_newestfile(self.outPath), self.info)
            pair=None
            if Mem.MolecularSystem.instruction.link.__len__()==2:
                pair_1,cloest_1 = Mem.findPair(Mem.MolecularSystem.instruction.link[0])
                pair_2,cloest_2 = Mem.findPair(Mem.MolecularSystem.instruction.link[1])
                
                if (pair_1 and pair_2 and (cloest_1>cloest_2)) or ((not pair_1) and pair_2):
                    pair=pair_2
                    link=Mem.MolecularSystem.instruction.link[1]
                    
                elif (pair_1 and pair_2 and (cloest_1<cloest_2)) or ((not pair_2) and pair_1):
                    pair=pair_1
                    link=Mem.MolecularSystem.instruction.link[0]
                else:
                    pair=None
            elif Mem.MolecularSystem.instruction.link.__len__()==1:
                link=Mem.MolecularSystem.instruction.link[0]
                pair,_=Mem.findPair(link)
            if pair:
                self.tryNum = 1
                Mem.Update(pair,link=link)
                self.bondsNum += 1
                Mem.MolecularSystem.lmp.writeNewLmpData(
                    f'{self.outPath}/bond_{self.bondsNum}.lmps')
                return True
            else:
                print("[INFO]\tNot find Bond, next try.")
                # number of attempts
                if self.tryNum > self.trylimit:
                    print("Not Find Bond, Break, code=10")
                    return False
                if Mem.MolecularSystem.instruction.equilibrium:
                    self.get_newestfile(self.outPath, self.lmprundir)
                    if self.bondsNum < 1: # Bond count less than 1, cannot get membrane width, cannot perform equilibrium
                        self.run_lammps(self.minfile)
                        print("[INFO]\tNo Bond is created, Do not check Equilibrium ")
                        shutil.move(f'{self.lmprundir}/min.lmps',
                                    f'{self.outPath}/bond_{self.bondsNum}_{self.tryNum}.lmps')
                    else:
                        if Mem.MolecularSystem.equilibrium():
                            Mem.MolecularSystem.lmp.writeNewLmpData(f'{self.outPath}/bond_{self.bondsNum}_{self.tryNum}.lmps')
                            self.get_newestfile(self.outPath, self.lmprundir)
                            print("[INFO]\tEquilibrium MIN")
                            self.run_lammps(self.minfile)
                            shutil.move(f'{self.lmprundir}/min.lmps',
                                        f'{self.outPath}/bond_{self.bondsNum}_{self.tryNum}.lmps')
                            shutil.move(f'{self.lmprundir}/log.lammps',f'trajectory/bond_{self.bondsNum}_{self.tryNum}_log.lammps')
                
                self.get_newestfile(self.outPath, self.lmprundir)
                self.run_lammps(self.mdfile)
                shutil.move(f'{self.lmprundir}/md.lmps',
                                f'{self.outPath}/bond_{self.bondsNum}_{self.tryNum}.lmps')
                shutil.move(f'{self.lmprundir}/hydration.lammpstrj',f'trajectory/bond_{self.bondsNum}_{self.tryNum}.lammpstrj')
                shutil.move(f'{self.lmprundir}/log.lammps',f'trajectory/bond_{self.bondsNum}_{self.tryNum}_log.lammps')
            self.tryNum += 1
            print(f'-'*10)

    def Prepare(self):
        Mem = InterfaceReaction(f'data.lmps', self.info,
                                f'{self.outPath}/bond_0.lmps')
        Mem.initCharges()
        self.get_newestfile(self.outPath, self.lmprundir)

    def remove_charges(self):
        Mem = InterfaceReaction(self.find_newestfile(
            self.outPath), self.info, f'bond_final.lmps')
        Mem.finalDeal()

    def run_lammps(self, infile):
        lmprun = f"{self.lammps} -i ../{infile}"
        os.chdir(self.lmprundir)
        if 'min' in infile:
            print(f'[INFO]\tRunning LAMMPS, Energy Minimization...')
        elif 'md' in infile:
            print(f'[INFO]\tRunning LAMMPS, Dynamics...')
        else:
            print(f'[INFO]\tRunning LAMMPS, {infile}...')
        
        FailedCounts=0
        
        while True:
            try:
                result = subprocess.run(
                    lmprun, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                # Check if LAMMPS ran successfully
                if result.returncode != 0:
                    current_time = datetime.now()
                    formatted_time = current_time.strftime('%Y-%m-%d %H:%M:%S')
                    print(f'-----{formatted_time}-----')
                    if result.returncode == 11:  # Segmentation fault
                        FailedCounts+=1
                        if FailedCounts>3:
                            print(f'[ERROR]\tToo many failures {FailedCounts}, Stop trying!')
                            raise subprocess.CalledProcessError(result.returncode, result.args, output=result.stdout, stderr=result.stderr)
                        else:
                            print(f"[ERROR]\tLAMMPS segmentation fault detected! Retrying...[{FailedCounts}]")
                    else:
                        print(f"[ERROR]\tLAMMPS run failed!")
                        raise subprocess.CalledProcessError(result.returncode, result.args, output=result.stdout, stderr=result.stderr)
                else:
                    break  # If run succeeded, exit loop

            except subprocess.CalledProcessError as e:
                print(f"[ERROR]\tCommand '{e.cmd}' returned non-zero exit status {e.returncode}.")
                print(f"[ERROR]\tStandard output:\n{e.output}")
                print(f"[ERROR]\tStandard error:\n{e.stderr}")
                if e.returncode != 11:
                    break  # If not segmentation fault, exit loop
            finally:
                os.chdir("..")

    def get_newestfile(self, folder_path: str, target_path: str) -> None:
        filestr = self.find_newestfile(folder_path)
        shutil.copy(f"{filestr}", f'{target_path}/data.lmps')

    def find_newestfile(self, folder_path: str) -> str:
        def version_key(s: str):
            if (not s.startswith("bond_")) and (not s.endswith(".lmps")):
                raise KeyError(
                    f"outdir have file which not a 'bond_*.lmps' ,{s}")
            s = list(map(int, s[5:-5].split('_')))
            if len(s) == 1:
                major = s[0]
                minor = -1
            else:
                major, minor = s
            return (major, minor)
        files = os.listdir(folder_path)
        max_file = max(files, key=version_key)
        return f"{folder_path}/{max_file}"
    
    def get_max_bondsNum_tryNum(self,outPath:str) -> tuple[int, int]:
        # Regular expression to match filename
        pattern = re.compile(r"bond_(\d+)(?:_(\d+))?\.lmps")
        
        bondsNum_tryNum_dict = {}
        
        # Iterate through all files in outPath directory
        for filename in os.listdir(outPath):
            match = pattern.match(filename)
            if match:
                # Extract bondsNum and tryNum
                bondsNum = int(match.group(1))
                tryNum = int(match.group(2)) if match.group(2) else 0
                
                # If bondsNum exists, update its maximum tryNum; otherwise, add to dictionary
                if bondsNum in bondsNum_tryNum_dict:
                    bondsNum_tryNum_dict[bondsNum] = max(bondsNum_tryNum_dict[bondsNum], tryNum)
                else:
                    bondsNum_tryNum_dict[bondsNum] = tryNum
        
        # If dictionary is not empty, get maximum bondsNum and corresponding maximum tryNum
        if bondsNum_tryNum_dict:
            max_bondsNum:int = max(bondsNum_tryNum_dict.keys())
            max_tryNum:int = bondsNum_tryNum_dict[max_bondsNum]
        else:
            max_bondsNum = 0
            max_tryNum = 0
        
        return max_bondsNum, max_tryNum

    def readconfig(self):

        with open('config/sample.yaml') as file:
            data = yaml.safe_load(file)

        self.bondsNum = data.get("bondsNum", 0)
        self.bondLimits = data.get("bondLimits", 600)
        self.tryNum = data.get("tryNum", 1)
        self.trylimit = data.get("trylimit", 100)
        self.outPath = data.get("outPath", "outs")
        self.minfile = data.get("minfile", "min.in")
        self.mdfile = data.get("mdfile", "md.in")
        self.lammps = data.get("lammps", None)
        self.lmprundir = data.get("lmprundir", "lammpsrun")

        self.info = {"typedata": readConfig.types(
            data), "instruction": readConfig.instruction(data)}
        
        # self.bondsNum, self.tryNum = self.get_max_bondsNum_tryNum(self.outPath)


if __name__ == '__main__':
    loop = membranepolymer()
