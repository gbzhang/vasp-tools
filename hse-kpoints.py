from pylab import *
import os

ibzkpt = "../1-pbe/IBZKPT"
os.system("cat " + ibzkpt + " > KPOINTS")

f = open("../2-pbe-bands/OUTCAR", 'r')

while True:
	line = f.readline()
	if "k-points in reciprocal lattice and weights" in line:
		nl = f.readline()
		while nl != " \n":
			l = " ".join(nl.split()[0:3]) + "  0.0"
			os.system("echo " + l + " >> KPOINTS")
			nl = f.readline()
		break
