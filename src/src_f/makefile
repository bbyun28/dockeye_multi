FFLAGS = -g -C -Wall -ffixed-line-length-132
MODULE = ../dockeyeM_energy.so
PYFILE = ../dockeyeM_c.py
SRC_F = src_f
# on mac, better use the f2py that comes with the PyMol
#F2PYHOME = ~/miniconda3/bin/f2py
F2PYHOME = /Applications/PyMOL_py37.app/Contents/bin/
#
# check fortran code for compile errors
code:
	gfortran $(FFLAGS) -c dockeyeM_energy.f
# use simple f2py make
main:
	$(F2PYHOME)/f2py  -c --f77flags="-g" dockeyeM_energy.f -m dockeyeM_energy
	#python -m numpy.f2py -c dockeyeM_energy.f -m dockeyeM_energy
	#ln -s -f $(SRC_F)/dockeyeM_f_v1.0.py $(PYFILE)
	#ln -s -f $(SRC_F)/dockeyeM_energy.cpython-37m-darwin.so $(MODULE)

# use  smarter method with in/out intent statements in *.pyf file- not working right now
# but is it really needed?
smart:
	#python -m numpy.f2py dockeyeM_energy.f -m dockeyeM_energy -h dockeyeM_energy.pyf
	python -m numpy.f2py -c dockeyeM_energy.pyf dockeyeM_energy.f
