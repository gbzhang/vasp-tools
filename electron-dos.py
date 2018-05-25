from pylab import *
rcParams.update({'font.size': 48, 'text.usetex': True})

f = open('DOSCAR')
ef = float(f.readlines()[5].split()[3])
f.close()

dos = genfromtxt('DOSCAR', skip_header=6)

el = -10
eh = 10

xh = max(dos[(dos[:,0] < eh) & (dos[:,0] > el),1])
xh *= 1.5

f=figure(figsize=(8,12))
plot(dos[:,1],dos[:,0]-ef)
axhline(0,color='k',linewidth=2)
axvline(0,color='k',linewidth=1)
ylim(el, eh)
xlim(0,xh)
xlabel('edos')
ylabel("Energy (eV)")
#savefig('edos_mote2.png', bbox_inches='tight', dpi=200)
tight_layout()
show()
