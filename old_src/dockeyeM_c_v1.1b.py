#############################################
#  Author: Kim Sharp
#  Date:   7/21/2017
#  
# branch off dockeye_c_v2.3.py  to dock multiple conformations
# of ligand simultaneously
# usage inside pymol command window: de('protein_target.pdb','ligand.pdb')
# target protein 1st, then ligand- ligand pdb file must have at least one model
# bracketed by MODEL, ENDMDL records
# TODO: 
# display any one of conformers- 
# compute energy for all conformers, but only generate dockeye object
# for lowest energy conformer
#
#############################################
import sys
import math
from pymol.callback import Callback
from pymol.cgo import *
from pymol import cmd
#from pdb_methods import *
import dockeyeM_energy
#import numpy as np

#=======================================
# include pdb_methods so code selfcontained
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
    def readfile(self,pdb_name):
        try:
            #pdb_file = open(pdb_name,"r") # python3
            pdb_file = file(pdb_name,"r")  # python2
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
        for entry in contents:
            if (entry[0:4] == 'ATOM') or (entry[0:6] == 'HETATM'):
                self.natom +=1
                xyz[0] = float(entry[30:38])
                xyz[1] = float(entry[38:46])
                xyz[2] = float(entry[46:54])
                self.coords.append([xyz[0],xyz[1],xyz[2]])
                self.bfact.append(float(entry[61:67]))
                atname = entry[11:17]
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                if(entry[12] != 'H'):
                    atype = entry[12:14]
                    #atype = ' ' + entry[12]
                else:
                    atype = ' ' + entry[12]
                if ( self.atom_rads.has_key(atype)):
                    self.radius.append(self.atom_rads[atype])
                else:
                    print ("atom radius not in dictionary", atype)
                    self.radius.append(0.0)
    def readligand(self,pdb_name):
        #
        # read multiple entries, delimited by MODEL/ENDMDL records
        #
        try:
            #pdb_file = open(pdb_name,"r") # python3
            pdb_file = file(pdb_name,"r")  # python2
            #print 'opening file'
        except:
            # print "can't find file", pdb_name
            sys.exit("can't find file")
        contents = pdb_file.readlines()
        pdb_file.close()
        #print(type(contents))
        print('lines read: ', len(contents))
        self.natom = 0
        self.nmodel = 0
        xyz = [0.,0.,0.]
        nline = 0
        while(nline < len(contents)):
          entry = contents[nline]
          while(entry[0:5] != 'MODEL'):
            nline += 1
            if(nline == len(contents)):break
            entry = contents[nline]
          if(nline == len(contents)):break
          #print('found model: ',nline)
          self.nmodel += 1
          #
          # if first model, create temporary pdb file for pymol
          if(self.nmodel == 1):
            tmp_pdb = open('dockeye_lig.pdb','w')
          nline += 1
          entry = contents[nline]
          print(entry)
          while(entry[0:6] != 'ENDMDL'):
            if (entry[0:4] == 'ATOM') or (entry[0:6] == 'HETATM'):
              #
              xyz[0] = float(entry[30:38])
              xyz[1] = float(entry[38:46])
              xyz[2] = float(entry[46:54])
              self.coords.append([xyz[0],xyz[1],xyz[2]])
              # only do names, charge, radius for first model
              if(self.nmodel == 1):
                tmp_pdb.write(entry)
                self.natom +=1
                self.bfact.append(float(entry[61:67]))
                atname = entry[11:17]
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                if(entry[12] != 'H'):
                  atype = entry[12:14]
                  #atype = ' ' + entry[12]
                else:
                  atype = ' ' + entry[12]
                if ( self.atom_rads.has_key(atype)):
                  self.radius.append(self.atom_rads[atype])
                else:
                  print ("atom radius not in dictionary", atype)
                  self.radius.append(0.0)
              # end if model 1
            # end if atom entry
            nline += 1
            entry = contents[nline]
          # end of model
          if(self.nmodel == 1):
            tmp_pdb.close()
        # end of file
        if(self.nmodel == 0):
          print('ERROR: ligand file must have at least one conformer')
          print('bracketed by MODEL, ENDMDL records')
          sys.exit()
        else:
          print('# of models: ',self.nmodel)
          nrec = len(self.coords)
          nrec_need = self.nmodel*self.natom
          if(nrec == nrec_need):
            print('found all coordinates: ',nrec_need)
          else:
            print('mismatch between expected & found # of coords ',nrec_need,nrec)
            sys.exit()
#=======================================
# defs to create cgo object directly rather than reading from file
def pymol_cgo_new(ctype):
    cgo_obj = [LINEWIDTH, 3.0,BEGIN,ctype]
    return cgo_obj
def pymol_cgo_end(cgo_obj):
    cgo_obj.append(END)
def pymol_cgo_addline(cgo_obj,vbeg,vend,cbeg,cend):
    #print 'in addline: ',vbeg,vend
    # start color
    rbeg = 1.
    gbeg = 1.
    bbeg = 1.
    if((cbeg >= 0.) and (cbeg < 0.5)):
      rbeg = 0.
      gbeg = 2.*cbeg
      bbeg = (1. - 2.*cbeg)
    elif((cbeg >= 0.5) and (cbeg <= 1.0)):
      rbeg = 2.*cbeg  - 1.
      gbeg = 2. - 2.*cbeg
      bbeg = 0.
    elif((cbeg > 1.0) and (cbeg <= 1.5)):
      rbeg = 1.
      gbeg = 0.
      bbeg = 2.*(cbeg - 1.0)
    # end color
    rend = 1.
    gend = 1.
    bend = 1.
    if((cend >= 0.) and (cend < 0.5)):
      rend = 0.
      gend = 2.*cend
      bend = (1. - 2.*cend)
    elif((cend >= 0.5) and (cend <= 1.0)):
      rend = 2.*cend  - 1.
      gend = 2. - 2.*cend
      bend = 0.
    elif((cend > 1.0) and (cend <= 1.5)):
      rend = 1.
      gend = 0.
      bend = 2.*(cend - 1.0)
    temp = [COLOR,rbeg,gbeg,bbeg,VERTEX,vbeg[0],vbeg[1],vbeg[2],
            COLOR,rend,gend,bend,VERTEX,vend[0],vend[1],vend[2]]
    #print 'temp', temp
    for i  in range(len(temp)):
      cgo_obj.append(temp[i])
def pymol_cgo_addtri(cgo_obj,v1,v2,v3,color,nm):
    red = 1.
    green = 1.
    blue = 1.
    coords = [v1,v2,v3]
    temp = []
    for k in range(3):
      if((color[k] >= 0.) and (color[k] < 0.5)):
        red = 0.
        green = 2.*color[k]
        blue = (1. - 2.*color[k])
      elif((color[k] >= 0.5) and (color[k] <= 1.0)):
        red = 2.*color[k]  - 1.
        green = 2. - 2.*color[k]
        blue = 0.
      elif((color[k] > 1.0) and (color[k] <= 1.5)):
        red = 1.
        green = 0.
        blue = 2.*(color[k] - 1.0)
      temp1 = [NORMAL,nm[0],nm[1],nm[2],COLOR,red,green,blue,VERTEX,coords[k][0],coords[k][1],coords[k][2]]
      for i in range(len(temp1)):
        temp.append(temp1[i])
    #print 'temp', temp
    for i  in range(len(temp)):
      cgo_obj.append(temp[i])
#=======================================
# this is the def executed on line in pymol
#def de(pdbfile1="ag.pdb",pdbfile2="ab.pdb",charges=True):
def de(pdbfile1="ab.pdb",pdbfile2="ag.pdb",charges=True,logscale=True,dielectric=80.,eps=0.1):
#def de(pdbfile1="ag.pdb",pdbfile2="ligand.pdbqt",charges=True,logscale=True,dielectric=80.,eps=0.1):
  # extract names, and create the 'object' name in the pymol menu window
  pdbobj1 = pdbfile1[:-4]
  cmd.load(pdbfile1, pdbobj1)
  #pdbobj2 = pdbfile2[:-4]
  #cmd.load(pdbfile2, pdbobj2)
  #
  # Dockeye class reads pdb file upon init
  pdbobj2 = 'dockeye_lig'
  obj = Dockeye(pdbfile1,pdbfile2,pdbobj1,pdbobj2,charges,logscale,dielectric,eps)
  #
  # but pymol also loads first ligand model for display
  cmd.load('dockeye_lig.pdb',pdbobj2)
  os.system('rm -f dockeye_lig.pdb')
  cmd.zoom()
  cmd.do("set auto_zoom, off",echo=0)
  name = 'dockeye'
  cmd.load_callback(obj,name)
  return obj

#=======================================
class Dockeye(Callback):
  def __init__(self,pdbfile1,pdbfile2,pdbobj1,pdbobj2,charges,logscale,dielectric,eps):
    # calling arguments
    self.pdbfile1 = pdbfile1
    self.pdbfile2 = pdbfile2
    self.pdbobj1 = pdbobj1
    self.pdbobj2 = pdbobj2
    self.logscale = logscale
    self.dielectric = dielectric
    self.eps = eps
    print('Initializing Dockeye...')
    print('pdb file 1: ',pdbfile1)
    print('pdb file 2: ',pdbfile2)
    print('charges, logscale energy bars: ',charges,logscale)
    print('energy parameters dielectric: %8.3f VDW depth: %8.3f' % (dielectric,eps))
    # read original pdb coords
    self.pdb1 = pdb_struct()
    self.pdb1.readfile(self.pdbfile1)
    self.pdb2 = pdb_struct()
    self.pdb2.readligand(self.pdbfile2)
    self.iconf = self.pdb2.nmodel
    #self.iconf = 2
    self.objmat1 = [1.,0.,0.,1., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
    self.objmat2 = [1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
    self.my_view = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
    #print name
    #print self.pdb1.coords
    #print self.pdb2.coords
    self.gcen1 = [0.,0.,0.]
    qtot1 = 0.
    qtot2 = 0.
    self.energy = [0.,0.,0.]
    self.energy_min = 0.
    for i in range(self.pdb1.natom):
      for k in range(3):
        self.gcen1[k] += self.pdb1.coords[i][k]
    for k in range(3):
      self.gcen1[k] /= self.pdb1.natom
    self.gcen2 = [0.,0.,0.]
    for i in range(self.pdb2.natom):
      for k in range(3):
        self.gcen2[k] += self.pdb2.coords[i][k]
    for k in range(3):
      self.gcen2[k] /= self.pdb2.natom
    if(not charges):
      for i in range(self.pdb1.natom):
        self.pdb1.bfact[i] = 0.
      for i in range(self.pdb2.natom):
        self.pdb2.bfact[i] = 0.
    if(charges):
      for i in range(self.pdb1.natom):
        qtot1 += self.pdb1.bfact[i]
      for i in range(self.pdb2.natom):
        qtot2 += self.pdb2.bfact[i]
    print '# of atoms 1: %6d   2: %6d' % (self.pdb1.natom,self.pdb2.natom)
    print 'geometric centers: '
    print '1:  %8.3f %8.3f %8.3f ' % (self.gcen1[0],self.gcen1[1],self.gcen1[2])
    print '2:  %8.3f %8.3f %8.3f ' % (self.gcen2[0],self.gcen2[1],self.gcen2[2])
    print 'net charge 1:  %8.3f 2:  %8.3f ' % (qtot1,qtot2)
    #
    # open timestamped logfile
    #
    os.system('date > dockeye_date.tmp')
    date_tmp = open('dockeye_date.tmp')
    date_raw = date_tmp.readline()
    date_tmp.close()
    os.system('rm -f dockeye_date.tmp')
    #print(date_raw)
    df = date_raw.split()
    #print(df)
    file_name = 'dockeye_' + df[2] + df[1] + df[5] + '_' + df[3] + '.log'
    print('writing to logfile: ',file_name)
    self.dockeye_log = open(file_name,'w')
    self.dockeye_log.write('pdbfile 1: '+ pdbfile1+'\n')
    self.dockeye_log.write('pdbfile 2: '+ pdbfile2+'\n')
    self.dockeye_log.write( '# of atoms 1: %6d   2: %6d\n' % (self.pdb1.natom,self.pdb2.natom))
    self.dockeye_log.write( 'geometric centers: \n')
    self.dockeye_log.write( '1:  %8.3f %8.3f %8.3f \n' % (self.gcen1[0],self.gcen1[1],self.gcen1[2]))
    self.dockeye_log.write( '2:  %8.3f %8.3f %8.3f \n' % (self.gcen2[0],self.gcen2[1],self.gcen2[2]))
    self.dockeye_log.write( 'net charge 1:  %8.3f 2:  %8.3f \n' % (qtot1,qtot2))
    self.dockeye_log.write('energy parameters dielectric: %8.3f VDW depth: %8.3f\n' % (dielectric,eps))
    self.dockeye_log.write('# of ligand conformers: %6d\n' % (self.pdb2.nmodel))


  def __call__(self):
        # get view on screen
        my_view = cmd.get_view()
        delta_mv = 0.
        nbest = 1
        for i in range(18):
          delta_mv = max(delta_mv,abs(my_view[i] - self.my_view[i]))
          self.my_view[i] = my_view[i] 
        #print my_view
        # get orientation/position matrices for two molecules
        # how does pymol define rotation center of molecule? 
        # apparnetly by geometric average
        pdbmat1 = cmd.get_object_matrix(self.pdbobj1)
        pdbmat2 = cmd.get_object_matrix(self.pdbobj2)
        delta_mm = 0.
        for i in range(12):
          delta_mm = max(delta_mm,abs(pdbmat1[i] - self.objmat1[i]))
          delta_mm = max(delta_mm,abs(pdbmat2[i] - self.objmat2[i]))
          self.objmat1[i] = pdbmat1[i] 
          self.objmat2[i] = pdbmat2[i] 
        if(delta_mm > 0.01): # we only do expensive energy calc if pose changed
          do_mm = True
        else:
          do_mm = False
        if((delta_mv > 0.01) or do_mm): # we only update if pose or view changed
          cgo_obj = pdb_interaction(pdbmat1,pdbmat2,self.pdb1,self.pdb2,self.gcen1,self.gcen2,
                  self.energy,do_mm,self.logscale,self.dielectric,self.eps,nbest)
          #
          # write new best pose to logfile
          #
          et = self.energy[0]
          ee = self.energy[1]
          ev = self.energy[2]
          if(self.energy[0] < self.energy_min):
            print('             NEW MIN ee: %12.3g ev: %12.3g et: %12.3g' % (ee,ev,et))
            self.energy_min = et
            self.dockeye_log.write('new min: %12.5g %12.5g %12.5g \n' % (ee,ev,et))
            for i in range(4):
              for j in range(4):
                indx = j + 4*i
                self.dockeye_log.write('%12.5f ' % (pdbmat1[indx]))
              self.dockeye_log.write('\n')
            for i in range(4):
              for j in range(4):
                indx = j + 4*i
                self.dockeye_log.write('%12.5f ' % (pdbmat2[indx]))
              self.dockeye_log.write('\n')
          else:
            if(do_mm): 
              print('Current energy: ee: %12.3g ev: %12.3g et: %12.3g' % (ee,ev,et))
        if(do_mm):
          cmd.delete('dockeye_obj')
          cmd.load_cgo(cgo_obj,'dockeye_obj')
        if(self.iconf != 0):
          draw_ligand(pdbmat2,self.pdb2,self.gcen2,self.iconf)
          self.iconf = 0

#=======================================
def rot_vec(rmt,vec,inv=0):
  vec_rot = [0.,0.,0.]
  if(inv == 0):
    for i in range(3):
      for j in range(3):
        vec_rot[i] += rmt[i][j]*vec[j]
  else:
    for i in range(3):
      for j in range(3):
        vec_rot[i] += rmt[j][i]*vec[j]
  return vec_rot

def vdot(v1,v2):
  dot = 0.
  for k in range(3):
    dot += v1[k]*v2[k]
  return dot
#
def vcross(v1,v2):
  cross = [0.,0.,0.]
  cross[0] = v1[1]*v2[2] - v1[2]*v2[1]
  cross[1] = v1[2]*v2[0] - v1[0]*v2[2]
  cross[2] = v1[0]*v2[1] - v1[1]*v2[0]
  return cross
#
def vnorm(v1):
  dot = vdot(v1,v1)
  dot = math.sqrt(dot)
  if(dot > 0.):
    for k in range(3):
      v1[k] /= dot
  return dot
#
def vperp(v1):
  v2 = [v1[1],v1[2],v1[0]]
  v3 = vcross(v1,v2)
  return v3

#=======================================
def pnl_make(rmt1,rmt2,gcen1,gcen2,trn1,trn2,energy):
  """
  refresh and display panel with energy bars and orientations
  """
  my_view = cmd.get_view()
  mod_center = [my_view[12],my_view[13],my_view[14]]
  rmtView = [[my_view[0], my_view[3], my_view[6]],
    [my_view[1], my_view[4], my_view[7]],
    [my_view[2], my_view[5], my_view[8]]] 
  cam_center = [my_view[9],my_view[10],my_view[11]]
  d_front = my_view[15]
  d_back = my_view[16]
  scale = d_back - d_front
  xsize = abs(cam_center[2])/22.
  xmove = 3.0
  #
  # create pnl objects
  #
  #
  # energy bars
  # +ve energies on log scale, -ve on linear scale now
  #
  ltype = TRIANGLES
  bar_obj = pymol_cgo_new(ltype)
  dy = -0.2
  #ecut = 0.5
  ecut = 0.2
  et = energy[0]
  ee = energy[1]
  ev = energy[2]
  if(et < -ecut):
    et_color = 0. # blue
    #et_size = -1.*xsize
    # log scale
    #et_size = (0.1 -1.*math.log10(-et/ecut))*xsize*2.
    # linear scale
    et_size = 0.5*et*xsize
  elif(et > ecut):
    et_color = 1. # red
    #et_size = 1.*xsize
    et_size = (0.1 + 1.*math.log10(et/ecut))*xsize*2.
  else:
    et_color = 2. # white
    et_size = 0.1*xsize
  #
  if(ee < -ecut):
    ee_color = 0. # blue
    #ee_size = -1.*xsize
    #ee_size = (0.1 -1.*math.log10(-ee/ecut))*xsize*2.
    # linear scale
    ee_size = 0.5*ee*xsize
  elif(ee > ecut):
    ee_color = 1. # red
    #ee_size = 1.*xsize
    ee_size = (0.1 + 1.*math.log10(ee/ecut))*xsize*2.
  else:
    ee_color = 2. # white
    ee_size = 0.1*xsize
  #
  if(ev < -ecut):
    ev_color = 0. # blue
    #ev_size = -1.*xsize
    #ev_size = (0.1 -1.*math.log10(-ev/ecut))*xsize*2.
    # linear scale
    ev_size = 0.5*ev*xsize
  elif(ev > ecut):
    ev_color = 1. # red
    #ev_size = 1.*xsize
    ev_size = (0.1 + 1.*math.log10(ev/ecut))*xsize*2.
  else:
    ev_color = 2. # white
    ev_size = 0.1*xsize
  #
  # at left, vertical bars
  #
  dx = -0.2 # bar separation
  #-----------------
  # E total
  #-----------------
  dend = [0.85*dx*xsize,0.,0.]
  dend_rot = rot_vec(rmtView,dend,inv=1)
  nm = [0.,0.,1.]
  nm_rot = rot_vec(rmtView,nm,inv=1)
  end1 = [0.,0.,0.]
  beg1 = [0.,0.,0.]
  #
  et_beg = [0.,0.,0.]
  et_end = [0.,et_size,0.]
  et_end_rot = rot_vec(rmtView,et_end,inv=1)
  et_offset = [-1.3*xmove*xsize,0.,0.]
  et_offset_rot = rot_vec(rmtView,et_offset,inv=1)
  for k in range(3):
    et_beg[k] = et_offset_rot[k] + mod_center[k]
    et_end_rot[k] += et_offset_rot[k] + mod_center[k]
    end1[k] = et_end_rot[k] - dend_rot[k]
    beg1[k] = et_beg[k] - dend_rot[k]
  color = [2.,et_color,et_color]
  pymol_cgo_addtri(bar_obj,et_beg,et_end_rot,end1,color,nm_rot)
  color = [2.,2.,et_color]
  pymol_cgo_addtri(bar_obj,et_beg,beg1,end1,color,nm_rot)
  #-----------------
  # E electrostatic
  #-----------------
  ee_beg = [0.,0.,0.]
  ee_end = [0.,ee_size,0.]
  ee_end_rot = rot_vec(rmtView,ee_end,inv=1)
  ee_offset = [(-1.3*xmove+2.*dx)*xsize,0.,0.]
  ee_offset_rot = rot_vec(rmtView,ee_offset,inv=1)
  for k in range(3):
    ee_beg[k] = ee_offset_rot[k] + mod_center[k]
    ee_end_rot[k] += ee_offset_rot[k] + mod_center[k]
    end1[k] = ee_end_rot[k] - dend_rot[k]
    beg1[k] = ee_beg[k] - dend_rot[k]
  color = [2.,ee_color,ee_color]
  pymol_cgo_addtri(bar_obj,ee_beg,ee_end_rot,end1,color,nm_rot)
  color = [2.,2.,ee_color]
  pymol_cgo_addtri(bar_obj,ee_beg,beg1,end1,color,nm_rot)
  #-----------------
  # E vdw
  #-----------------
  ev_beg = [0.,0.,0.]
  ev_end = [0.,ev_size,0.]
  # apply view angle
  ev_end_rot = rot_vec(rmtView,ev_end,inv=1)
  ev_offset = [(-1.3*xmove+dx)*xsize,0.,0.]
  ev_offset_rot = rot_vec(rmtView,ev_offset,inv=1)
  for k in range(3):
    ev_beg[k] = ev_offset_rot[k] + mod_center[k]
    ev_end_rot[k] += ev_offset_rot[k] + mod_center[k]
    end1[k] = ev_end_rot[k] - dend_rot[k]
    end1[k] = ev_end_rot[k] - dend_rot[k]
    beg1[k] = ev_beg[k] - dend_rot[k]
  color = [2.,ev_color,ev_color]
  pymol_cgo_addtri(bar_obj,ev_beg,ev_end_rot,end1,color,nm_rot)
  color = [2.,2.,ev_color]
  pymol_cgo_addtri(bar_obj,ev_beg,beg1,end1,color,nm_rot)
  #
  # finish up & display energy bars
  #
  pymol_cgo_end(bar_obj)
  cmd.delete('bar_obj')
  cmd.load_cgo(bar_obj,'bar_obj')

#=======================================
def pdb_interaction(pdbmat1,pdbmat2,pdb1,pdb2,gcen1,gcen2,energy,do_mm,logscale,dielectric,eps,
                    nbest):
  """
  # this is where we calculated interaction energy between
  # two molecules
  #
  # extract rotation matrices and translation vectors
  # note that if center of molecule not at origin
  # then rotation of a molecule moves the origin relative to its
  # local center, and that this displacment appears in translation so
  # this apparent trans needs to be removed to get actual translation
  """
  #=========================================
  # global parameters, these are the nominal values used by C
  # subroutine, energies are then post-scaled for actual parameters given as arguments to de()
  DIEL = 80.
  efact = 332./DIEL # dielectric constant factor gives kcal/mole for p+ unit charge, Angstroms
  EPS = 0.1 # depth of vdw potl. kcal/mole
  #=========================================
  # extract rot mat
  rmt1 = [[pdbmat1[0],pdbmat1[1],pdbmat1[2]],
         [pdbmat1[4],pdbmat1[5],pdbmat1[6]],
         [pdbmat1[8],pdbmat1[9],pdbmat1[10]]]
  rmt2 = [[pdbmat2[0],pdbmat2[1],pdbmat2[2]],
         [pdbmat2[4],pdbmat2[5],pdbmat2[6]],
         [pdbmat2[8],pdbmat2[9],pdbmat2[10]]]
  #
  # extract trans
  trn1 = [pdbmat1[3],pdbmat1[7],pdbmat1[11]]
  trn2 = [pdbmat2[3],pdbmat2[7],pdbmat2[11]]
  # find rotated origin in molc. local coords
  gcen1_rot = rot_vec(rmt1,gcen1)
  gcen2_rot = rot_vec(rmt2,gcen2)
  # displacement of origin subtracted from apparant trans to get actual trans
  for k in range(3):
    trn1[k] = trn1[k] - (gcen1[k] - gcen1_rot[k])
    trn2[k] = trn2[k] - (gcen2[k] - gcen2_rot[k])
  #
  # refresh display panel
  if(not do_mm):
    pnl_make(rmt1,rmt2,gcen1,gcen2,trn1,trn2,energy)
    return
  #
  #
  xyz = [0.,0.,0.]
  atom_data = [0.,0.] # stuff # atoms (as floats), coords, radii and charges in one long array to pass to energy_c
  nat1 = pdb1.natom
  nat2 = pdb2.natom
  #
  # molecule 1
  #
  for i in range(nat1):
    # apply rotations and translations
    for k in range(3):
      xyz[k] = pdb1.coords[i][k] - gcen1[k]
    xyz1 = rot_vec(rmt1,xyz)
    for k in range(3):
      xyz[k] = xyz1[k] + gcen1[k] + trn1[k]
    # store new coords
    for k in range(3):
      atom_data.append(xyz[k])
    atom_data.append(pdb1.radius[i])
    atom_data.append(pdb1.bfact[i])
  #
  # molecule 2
  #
  for i in range(nat2):
    # apply rotations and translations
    for k in range(3):
      xyz[k] = pdb2.coords[i][k] - gcen2[k]
    xyz2 = rot_vec(rmt2,xyz)
    for k in range(3):
      xyz[k] = xyz2[k] + gcen2[k] + trn2[k]
    # store new coords
    for k in range(3):
      atom_data.append(xyz[k])
    atom_data.append(pdb2.radius[i])
    atom_data.append(pdb2.bfact[i])
  #
  #==========================================
  # C subroutine version
  #==========================================
  atom_data[0] = float(nat1) # now we know # of atoms, put in front of data array
  atom_data[1] = float(nat2)
  energy_obj = dockeyeM_energy.energy_c(atom_data)
  ndata = int(energy_obj[0])
  # slice energy terms off end of data
  #energy[0] = energy_obj[ndata-3]
  energy[1] = energy_obj[ndata-2]*DIEL/dielectric
  energy[2] = energy_obj[ndata-1]*eps/EPS
  energy[0] = energy[1] + energy[2]
  energy_obj[0] = LINEWIDTH # after we extract length of data, put real cgo 1st keyword back
  # generate true Pymol object- one returned from C doesn't seem to work
  cgo_obj = []
  for i in range(ndata):
    cgo_obj.append(energy_obj[i])
  pnl_make(rmt1,rmt2,gcen1,gcen2,trn1,trn2,energy)
  return cgo_obj

def draw_ligand(pdbmat2,pdb2,gcen2,iconf):
  """
   if best ligand coformer has changed, regenerate the object
  """
  if(iconf == 0):
    return
  #=============================
  ibeg = (iconf - 1)*pdb2.natom
  ifin = ibeg + pdb2.natom
  ltype = LINES
  ligand_obj = pymol_cgo_new(ltype)
  print('drawing conf: ',iconf,ibeg,ifin,pdb2.natom)
  #
  vbeg = [0.,0.,0.]
  vend = [0.,0.,0.]
  cbeg = .5
  cend = .9
  #
  # extract rot mat
  rmt2 = [[pdbmat2[0],pdbmat2[1],pdbmat2[2]],
         [pdbmat2[4],pdbmat2[5],pdbmat2[6]],
         [pdbmat2[8],pdbmat2[9],pdbmat2[10]]]
  #
  # extract trans
  trn2 = [pdbmat2[3],pdbmat2[7],pdbmat2[11]]
  #
  # find rotated origin in molc. local coords
  gcen2_rot = rot_vec(rmt2,gcen2)
  #
  # displacement of origin subtracted from apparant trans to get actual trans
  for k in range(3):
    trn2[k] = trn2[k] - (gcen2[k] - gcen2_rot[k])
  #
  # extract and transform coords from required model
  xyz = [0.,0.,0.]
  crd2 = []
  for i in range(ibeg,ifin):
    # apply rotations and translations
    for k in range(3):
      xyz[k] = pdb2.coords[i][k] 
    for k in range(3):
      xyz[k] = pdb2.coords[i][k] - gcen2[k]
    xyz2 = rot_vec(rmt2,xyz)
    for k in range(3):
      xyz[k] = xyz2[k] + gcen2[k] + trn2[k]
    crd2.append([xyz[0],xyz[1],xyz[2]])
  #
  # generate bonds
  nbond = 0
  margin = 1.4
  for i in range(pdb2.natom):
    if(pdb2.name[i][2:3] == 'H'):
      radi = 1.1
    else:
      radi = pdb2.radius[i]
    for j in range(i+1,pdb2.natom):
#        print 'checking: ',i,j
        dist2 = 0.
        if(pdb2.name[j][2:3] == 'H'):
          radj = 1.1
        else:
          radj = pdb2.radius[j]
        for k in range(3):
            dist2 += (crd2[i][k] - crd2[j][k])**2
        dist = math.sqrt(dist2)
        overlap = dist + margin - radi - radj
        if (overlap < 0.):
            nbond +=1
            for k in range(3):
              vbeg[k] = crd2[i][k]
              vend[k] = crd2[j][k]
            if(pdb2.bfact[i] > 0.1):
              cbeg = 0.1
            elif(pdb2.bfact[i] < -0.1):
              cbeg = 0.9
            else:
              cbeg = 2.
            if(pdb2.bfact[j] > 0.1):
              cend = 0.1
            elif(pdb2.bfact[j] < -0.1):
              cend = 0.9
            else:
              cend = 2.
            pymol_cgo_addline(ligand_obj,vbeg,vend,cbeg,cend)
  print "number of bonds: ",nbond
  pymol_cgo_end(ligand_obj)
  cmd.delete('ligand_obj')
  cmd.load_cgo(ligand_obj,'ligand_obj')
