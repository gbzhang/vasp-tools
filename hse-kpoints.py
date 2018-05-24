from pylab import *
import os

ibzkpt = "../4-pbe-coarse/IBZKPT"
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
f.close()

f = open('KPOINTS', 'r')
lines = f.readlines()
f.close()


slines = []
[slines.append(line) for line in lines if len(line.split()) == 4 ]
nk = len(slines)
lines[1] = '     ' + str(nk) + '\n'

f = open('KPOINTS', 'w')
f.write("".join(lines))
f.close()

