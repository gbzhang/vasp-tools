import os
import sys
from pylab import *

def readPOSCAR(filedir):
	f = open(filedir + "/POSCAR", 'r')
	lines = f.readlines()
	# lines 1,2 are ignorable
	a = array([float(i) for i in lines[2].split()])
	b = array([float(i) for i in lines[3].split()])
	c = array([float(i) for i in lines[4].split()])

	species = lines[5].split()
	species_and_atomnums = lines[5] + lines[6]
	atomnums = [int(i) for i in lines[6].split()]
	natoms = sum(atomnums)
	
	positions = zeros([natoms, 3])

	# Line 7 should be "Direct"
	if lines[7].strip()[0] == 'D' or lines[7].strip()[0] == 'd':
		for i in range(natoms):
			positions[i] = array([float(j) for j in lines[8+i].split()])

	else:
		print("Cannot handle anything but Direct coordinates for now")
		exit(1)

	lat = array([a,b,c])
	
	return lat, positions, natoms, species_and_atomnums
	
def getInterpolatedPOSCAR(Ncurrent, Nmax, lat_i, lat_f, pos_i, pos_f, header_string, spec_atnum_i):
	poscar_string = header + "\n 1.00000\n"

	for j in range(3):
		lat_int = array(lat_i) + float(Ncurrent/Nmax)*(array(lat_f)-array(lat_i))
		poscar_string += " ".join(["%20.15f"%k for k in lat_int[j]]) + "\n"

	#print(poscar_string)
	poscar_string += spec_atnum_i 
	poscar_string += "Direct\n"
	positions = pos_i + float(Ncurrent/Nmax)*(pos_f - pos_i)
	natoms = shape(positions)[0]
	
	for j in range(natoms):
		poscar_string += " ".join(["%20.15f"%k for k in positions[j]]) + "\n"
	
	return poscar_string


def writePOSCARfiles(Nmax, lat_i, lat_f, pos_i, pos_f, header, spec_atnum_i):
	for n in range(1,Nmax):
		f = open(folders[n] + '/POSCAR', 'w')
		f.write(getInterpolatedPOSCAR(n, Nmax, lat_i, lat_f, pos_i, pos_f, header, spec_atnum_i))
		f.close()

	return 0

if __name__ == '__main__':
	header = 'SnSe' # Header for POSCAR files (must be the same for each)

	folders = sort([i for i in os.listdir('.') if '0' in i])
	Nmax = int(folders[-1])

	lat_i, pos_i, nat_i, spec_atnum_i = readPOSCAR(folders[0])
	lat_f, pos_f, nat_f, spec_atnum_f = readPOSCAR(folders[-1])

	if norm(lat_i - lat_f) > 1e-12:
		print("Warning: interpolating between two sets of lattice constants")

	writePOSCARfiles(Nmax, lat_i, lat_f, pos_i, pos_f, header, spec_atnum_i)
