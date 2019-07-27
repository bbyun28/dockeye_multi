# eventually figure out correct lib-dynload directory on each machine to insert link to *.so
#DYNLOAD = /home/kim/Downloads/pymol/ext/lib/python2.4/lib-dynload

# rules
.f.o:
	$(FCOMP) -c $(FFLAGS) $*.f

# home computer
home:
	gcc -O -fPIC -m32 -c dockeyeM_energy.c -o dockeyeM_energy.o -I/usr/include -I/home/kim/anaconda2/include/python2.7/
	gcc -m32 -shared -s dockeyeM_energy.o -lm  -L/home/kim/anaconda2/lib/python2.7/config/ -lpython2.7 -u MAIN__ -o dockeyeM_energy.so
# mac book air
macbook:
	gcc -g -shared -I/sw/include/python2.7 -lpython2.7 -L/sw/lib/python2.7/config/ -o dockeyeM_energy.so dockeyeM_energy.c
# office
office:
	gcc -O -fPIC -m32 -c dockeyeM_energy.c -o dockeyeM_energy.o -I/usr/include -I/home/sharp/anaconda2/include/python2.7/
	gcc -m32 -shared -s dockeyeM_energy.o -lm  -L/home/sharp/anaconda2/lib/python2.7/config/ -lpython2.7 -u MAIN__ -o dockeyeM_energy.so
# redtrooper
redtrooper:
	gcc -O -fPIC -c dockeyeM_energy.c -o dockeyeM_energy.o -I/usr/include -I/usr/include/python2.7/
	gcc -shared -s dockeyeM_energy.o -lm  -Lpython2.7 -lpython2.7 -u MAIN__ -o dockeyeM_energy.so
