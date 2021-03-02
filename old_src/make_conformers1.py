"""
read an auto dock format pdb file (*.pdbqt)
and generate conformers, write out in MODEL/ENDMDL format
for pymol and dockeye multi
"""
import sys
import math
#
# include pdb_methods so code selfcontained
# and we need special code for parsing the *.pdbqt ROOT/BRANCH
# torsion tree records
"""implement class and methods for handling pdb file"""
#import string
class pdb_struct:
    def __init__(self):
        """create some placeholders for coords, radius, atom name, atom type etc"""
        self.coords = []
        self.radius = []
        self.name = []
        self.res = []
        self.element = []
        self.resnum = []
        self.chain = []
        self.bfact = []
        self.atom_rads = {' C':1.8,' S':1.9,' O':1.6,' N':1.4,' P':1.8,' H':0.0,'ZN':1.4,
        ' Z':1.4,' B':2.46}
        self.root = [1,0]
    def readfile(self,pdb_name):
        try:
            pdb_file = open(pdb_name,"r") # python3
            #pdb_file = file(pdb_name,"r")  # python2
            #print 'opening file'
        except:
            # print "can't find file", pdb_name
            sys.exit("can't find file")
        contents = pdb_file.readlines()
        pdb_file.close()
        #print(type(contents))
        #print('lines read: ', len(contents))
        self.natom = 0
        xyz = [float(0.),float(0.),float(0.)]
        inroot = 0
        branch_stack = []
        nbranch = 0
        self.branch = []
        branch_atoms = []
        self.ntor = 0
        atom_name_number = {}
        self.tor_atoms = []
        for entry in contents:
            if(entry[0:7] == 'TORSDOF'):
              fields = entry.split()
              ntor_expt = int(fields[1])
              print ('expected # of torsions: ',ntor_expt)
              if((ntor_expt != nbranch) or (ntor_expt != self.ntor)):
                print('ERROR: found ',nbranch)
                sys.exit()
            if(entry[0:6] == 'REMARK'):
              if('between atoms' in entry):
                self.ntor += 1
                fields = entry.split()
                n1 = fields[5].split('_')
                name1 = n1[0]
                n2 = fields[7].split('_')
                name2 = n2[0]
                #print(name1,name2)
                self.tor_atoms.append([name1,name2])
            if(entry[0:4] == 'ROOT'):
              inroot = 1
            if(entry[0:7] == 'ENDROOT'):
              self.root[1] = self.natom
              print('found root: ',self.root)
              inroot = 1
            if(entry[0:6] == 'BRANCH'):
              nbranch += 1
              fields = entry.split()
              m1 = int(fields[1])
              m2 = int(fields[2])
              branch_atoms.append([m1,m2])
              branch_stack.insert(0,nbranch)
              self.branch.append([self.natom,0])
            if(entry[0:9] == 'ENDBRANCH'):
              i = branch_stack.pop(0) - 1
              self.branch[i][1] = self.natom
            if (entry[0:4] == 'ATOM') or (entry[0:6] == 'HETATM'):
                self.natom +=1
                xyz[0] = float(entry[30:38])
                xyz[1] = float(entry[38:46])
                xyz[2] = float(entry[46:54])
                self.coords.append([xyz[0],xyz[1],xyz[2]])
                self.bfact.append(float(entry[61:67]))
                atname = entry[11:17].strip()
                atom_name_number[atname] = self.natom
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                if(entry[12] != 'H'):
                    atype = entry[12:14]
                    #atype = ' ' + entry[12]
                else:
                    atype = ' ' + entry[12]
                #if ( self.atom_rads.has_key(atype)):
                if ( atype in self.atom_rads):
                    self.radius.append(self.atom_rads[atype])
                else:
                    print ("atom radius not in dictionary", atype)
                    self.radius.append(0.0)
    
        print('# of branches: ',nbranch)
        if(len(branch_stack) !=0):
          print('ERROR: branches: ',branch_stack)
        print('branches: ',self.branch)
        print('tor_atoms: ',self.tor_atoms)
        print('atom_name_number: ',atom_name_number)
        tor_numbers = []
        for i in range(self.ntor):
           j = atom_name_number[self.tor_atoms[i][0]]
           k = atom_name_number[self.tor_atoms[i][1]]
           if(j < k):
             tor_numbers.append([j,k])
           else:
             tor_numbers.append([k,j])
        print('tor_numbers : ',tor_numbers)
        print('branch atoms: ',branch_atoms)
        for i in range(self.ntor):
          a1 = branch_atoms[i][0]
          a2 = branch_atoms[i][1]
          print(i,'twisting about ',a1,a2)
          for j in range(a2,self.branch[i][1]):
            print('%4d ' % (j+1),end='')
          print(' \n')
# =======================
pdbfile = 'temp_model1.pdbqt'
#pdbfile = 'ag.pdbqt'
pdb1 = pdb_struct()
pdb1.readfile(pdbfile)
#
# now trace torsion tree
#
"""
for i in range(ntor):
  a1 = branch_atoms[i][0]
  a2 = branch_atoms[i][1]
  print(i,'twisting about ',a1,a2)
  for j in range(a2+1,branches[i][2]):
    print('%d4 ' % (j))
  print(' \n')
"""
