#############################################
#  Author: Kim Sharp
#  Date:   7/21/2018
# dockeyeM_getpose.py 
# extract pose from dockeye logfile
# add generation of movie file of ligand
# branch off dockeye_getpose.py to read models from input pdb file
# and best model out of logfile
#############################################
import sys
import math
from dockeye_methods import *
#=======================================
# MAIN
#=======================================
print('\n extract pose from dockeye logfile ')
print(' also generate movie file of ligand poses \n')
if(len(sys.argv) < 2):
  print('Usage: python dockeyeM_getpose.py dockeyeM_logfile')
  sys.exit()
file_name = sys.argv[1]
dockeye_log = open(file_name,'r')
contents = dockeye_log.readlines()
nline = len(contents)
dockeye_log.close()
print('lines read: ',nline)
#
# pdb files
#
iline = 0
fields = contents[iline].split()
pdbfile1 = fields[2]
iline += 1
fields = contents[iline].split()
pdbfile2 = fields[2]
print('pdb files used by dockeye: ')
print(pdbfile1,pdbfile2)
#
# get # of atoms
#
iline += 1
fields = contents[iline].split()
nat1 = int(fields[4])
nat2 = int(fields[6])
if(nat1 < nat2):
  print('warning: protein needs to be first pdb in file')
  sys.exit()
print('# atoms in protein, ligand: ',nat1,nat2)
#
# geometric centers
#
iline += 2
gcen1 = [0.,0.,0.]
fields = contents[iline].split()
for k in range(3):
  gcen1[k] = float(fields[k+1])
iline += 1
gcen2 = [0.,0.,0.]
fields = contents[iline].split()
for k in range(3):
  gcen2[k] = float(fields[k+1])
#print(gcen1,gcen2)
header = True
iline += 1
while(header):
  head = contents[iline][0:7]
  if(head == 'new min'):
    header = False
    istart = iline
    print('istart: ',istart)
  else:
    print(contents[iline][0:-1])
    iline += 1
#
#poses
#
npose = 0
et = []
ee = []
ev = []
nbest = []
while(iline < nline):
  head = contents[iline][0:7]
  #print(head)
  if(head == 'new min'):
    #print(head)
    fields = contents[iline].split()
    ee.append(float(fields[2]))
    ev.append(float(fields[3]))
    et.append(float(fields[4]))
    nbest.append(int(fields[6]))
    npose +=1
  iline += 1
print('# of poses: ',npose)
print(' ')
print('pose #      Ee           Ev           Et     model')
for i in range(npose):
  print('(%4d) %12.5f %12.5f %12.5f %6d' % (i+1,ee[i],ev[i],et[i],nbest[i]))
ipose = -1
xyz = [0.,0.,0.]
while((ipose <0) or(ipose > npose)):
  print('select pose # (1 - ',npose,')',)
  print(' enter 0 to exit')
  ipose = int(input('>> '))
  #ipose = npose # debug
if(ipose > 0):
  indx = 9*(ipose - 1) + istart  + 1
  imod = nbest[ipose - 1]
  print('pose, model: ',ipose,imod)
  pdbmat1 = []
  pdbmat2 = []
  for i in range(4):
    #print(contents[indx])
    fields = contents[indx].split()
    for k in range(4):
      pdbmat1.append(float(fields[k]))
    indx += 1
  for i in range(4):
    #print(contents[indx])
    fields = contents[indx].split()
    for k in range(4):
      pdbmat2.append(float(fields[k]))
    indx += 1
  #print(pdbmat1)
  #print(pdbmat2)
  #
  # read pdb files
  #
  pdb1 = pdb_struct()
  pdb1.readfile(pdbfile1)
  pdb2 = pdb_struct()
  pdb2.readligand(pdbfile2)
  print('# of atoms: ',pdb1.natom,pdb2.natom)
  #
  # transform
  #
  # extract rot mat
  rmt1 = [[pdbmat1[0],pdbmat1[1],pdbmat1[2]],
         [pdbmat1[4],pdbmat1[5],pdbmat1[6]],
         [pdbmat1[8],pdbmat1[9],pdbmat1[10]]]
  rmt2 = [[pdbmat2[0],pdbmat2[1],pdbmat2[2]],
         [pdbmat2[4],pdbmat2[5],pdbmat2[6]],
         [pdbmat2[8],pdbmat2[9],pdbmat2[10]]]
  #print(rmt1)
  #print(rmt2)
  #
  # extract trans
  trn1 = [pdbmat1[3],pdbmat1[7],pdbmat1[11]]
  trn2 = [pdbmat2[3],pdbmat2[7],pdbmat2[11]]
  #print(trn1)
  #print(trn2)
  #
  # find rotated origin in molc. local coords
  gcen1_rot = rot_vec(rmt1,gcen1)
  gcen2_rot = rot_vec(rmt2,gcen2)
  #
  # displacement of origin subtracted from apparant trans to get actual trans
  for k in range(3):
    trn1[k] = trn1[k] - (gcen1[k] - gcen1_rot[k])
    trn2[k] = trn2[k] - (gcen2[k] - gcen2_rot[k])
  #
  #
  # molecule 1
  #
  for i in range(pdb1.natom):
    # apply rotations and translations
    for k in range(3):
      xyz[k] = pdb1.coords[i][k] - gcen1[k]
    xyz1 = rot_vec(rmt1,xyz)
    for k in range(3):
      xyz[k] = xyz1[k] + gcen1[k] + trn1[k]
    for k in range(3):
      pdb1.coords[i][k] = xyz[k]
  #
  # molecule 2
  #
  for i in range(pdb2.natom):
    # apply rotations and translations
    i1 = i + imod*pdb2.natom
    for k in range(3):
      xyz[k] = pdb2.coords[i1][k] - gcen2[k]
    xyz2 = rot_vec(rmt2,xyz)
    for k in range(3):
      xyz[k] = xyz2[k] + gcen2[k] + trn2[k]
    for k in range(3):
      pdb2.coords[i][k] = xyz[k]
  #
  # outout pose pdb files
  #
  #
  # open and write atm file
  #
  if(ipose <10):
   numb = ('%1d'%(ipose))
  elif(ipose <100):
   numb = ('%2d'%(ipose))
  else:
   numb = ('%3d'%(ipose))
  #
  # first pdb
  #
  atm_file = pdbfile1[0:-4]+'_pose_'+numb+'.pdb'
  print('writing pose to ',atm_file)
  pdb_out = open(atm_file,'w')
  for i in range(pdb1.natom):
    #print('q: %8.3f  r: %8.3f ' % (pdb_charge[i],pdb.radius[i]))
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
    pdb1.chain[i], pdb1.resnum[i], pdb1.coords[i][0], \
    pdb1.coords[i][1],pdb1.coords[i][2],pdb1.radius[i],pdb1.bfact[i])
    pdb_out.write(string)
  pdb_out.close()
  #
  # second pdb
  #
  atm_file = pdbfile2[0:-4]+'_pose_'+numb+'.pdb'
  print('writing pose to ',atm_file)
  pdb_out = open(atm_file,'w')
  for i in range(pdb2.natom):
    #print('q: %8.3f  r: %8.3f ' % (pdb_charge[i],pdb.radius[i]))
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb2.name[i],pdb2.res[i], \
    pdb2.chain[i], pdb2.resnum[i], pdb2.coords[i][0], \
    pdb2.coords[i][1],pdb2.coords[i][2],pdb2.radius[i],pdb2.bfact[i])
    pdb_out.write(string)
  pdb_out.close()
#
# generate movie file of ligand poses
#
print('generating movie of poses...')
pdbl = pdb_struct()
gcenl = [0.,0.,0.]
pdbl.readligand(pdbfile2)
atm_file = pdbfile2[0:-4]+'_movie.pdb'
for k in range(3):
  gcenl[k] = gcen2[k]
pdb_out = open(atm_file,'w')
print('# of atoms: ',pdbl.natom)
iline = istart
nframe = 0
rmtl = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
trnl = [0.,0.,0.]
while(iline < nline):
  head = contents[iline][0:7]
  #print(head)
  if(head == 'new min'):
    fields = contents[iline].split()
    imod = int(fields[6])
    #print(head)
    nframe +=1
    #print ('f: ',nframe)
    pdb_out.write('MODEL %3d\n' % (nframe))
    iline += 4
    #
    # extract rot, trans for this frame
    #
    for i in range(3):
      iline +=1
      fields = contents[iline].split()
      trnl[i] = float(fields[3])
      for j in range(3):
        rmtl[i][j] = float(fields[j])
    #print(rmtl)
    #print(trnl)
    iline += 1
    #
    # find rotated origin in molc. local coords
    # then displacement of origin subtracted from apparant trans to get actual trans
    #
    gcenl_rot = rot_vec(rmtl,gcenl)
    for k in range(3):
      trnl[k] = trnl[k] - (gcenl[k] - gcenl_rot[k])
    for i in range(pdbl.natom):
      #
      # apply rotations and translations
      i1 = i + imod*pdbl.natom
      for k in range(3):
        xyz[k] = pdbl.coords[i1][k] - gcenl[k]
      xyz1 = rot_vec(rmtl,xyz)
      for k in range(3):
        xyz[k] = xyz1[k] + gcenl[k] + trnl[k]
      string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdbl.name[i],pdbl.res[i], \
      pdbl.chain[i], pdbl.resnum[i], xyz[0],xyz[1],xyz[2], \
      pdbl.radius[i],pdbl.bfact[i])
      pdb_out.write(string)
    # end of pdb write
  # next frame
  pdb_out.write('ENDMDL\n')
  iline += 1
print('# of movie frames written: ',nframe)
print(' to file: ',atm_file)
#pdb_out.close()
