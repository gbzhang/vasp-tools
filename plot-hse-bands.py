from pylab import *
rcParams.update({'font.size': 42, 'text.usetex': True})

def getEfermi():
	f = open("vasprun.xml",'r')
	while True:
		line = f.readline()
		if 'efermi' in line:
			ef = float(line.split()[2])
			break
	f.close()	
	return ef


def getNkpoints():
	f = open("../2-pbe-bands/EIGENVAL","r")
	for i in range(5):
		f.readline()
	nk = int(f.readline().split()[1])
	f.close()
	print("Found ", nk, " k points")
	return nk

def getbands(nk):
	f = open("EIGENVAL", "r")
	for i in range(5):
		f.readline()
	nbands = int(f.readline().split()[2])
	print("Found ", nbands, "bands")

	kp = zeros(nk)
	bands = zeros([nbands, nk])

	kvl = zeros(3)
	k = 0
	reached_zero_weight = False
	while True:
		line = f.readline()
		if line == '':
			break
		
		if '0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00' in line:
			reached_zero_weight = True

		if reached_zero_weight:
			line = line.split()
			if len(line) == 4:
				kv = array(list(map(float,line[0:3])))
				if k > 0:
					kp[k] = kp[k-1] + norm(kv-kvl)
				else:
					kp[k] = norm(kv)
				kvl = kv
				k += 1
			elif len(line) == 2:
				bands[0,k-1] = float(line[1])
				for i in range(1,nbands):
					bands[i, k-1] = float(f.readline().split()[1])

				f.readline() # blank
	f.close()

	# Get kpath from PBE KPOINTS file
	f = open("../2-pbe-bands/KPOINTS", "r")
	xtnames = []
	xtpos  = []
	k = 0
	for i in range(4):  #skip first 4 lines
		f.readline()
	while True:
		line = f.readline()
		if line == '':
			break
		else:
			xt = array(list(map(float,line.split()[0:3])))
			xtnames.append(line.split()[4])
			if k > 0:
				xtpos.append( xtpos[k-1] + norm(xt-xtl) )
			else:
				xtpos.append(norm(xt))
			k += 1
			xtl = xt
			f.readline()
			f.readline()
	f.close()
	xtpos.append(kp[-1])
	xtnames.append(xtnames[0])

	bnds = bands -ef
	gap = min(bnds[bnds>0.05])

	return kp, xtpos, xtnames, bnds, gap
	

if __name__ == '__main__':

	ef = getEfermi()
	nk = getNkpoints()
	
	kp, xtpos, xtnames, bands, gap = getbands(nk)
	
	figure()
	plot(kp, transpose(bands),color='k')
	axhline(0,color='g',linewidth=1)
	axhline(gap,color='r',linewidth=1)
	for i in range(1,len(xtpos)-1):
		axvline(xtpos[i],color='k',linewidth=1)
	ylabel('$E-E_f$ (eV)')
	xticks(xtpos,xtnames)
	ylim(-5,5)
	xlim(kp[0],kp[-1])
	title('gap=%f'%gap, fontsize=20)
	tight_layout()
	show()


