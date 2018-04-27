import yaml
from pylab import *
rcParams.update({'font.size': 42, 'text.usetex': True})

yam = yaml.load(open('band.yaml', 'r'))

# first get distance along k path
nb = len(yam['phonon'][0]['band'])  # number of bands
nq = yam['nqpoint']                 # number of q points

# values of q points
qp = array([yam['phonon'][i]['distance'] for i in range(nq)])

# values of frequencies in bands
freqs = zeros([nb, nq])
for j in range(yam['nqpoint']):
    freqs[:,j] = array([yam['phonon'][j]['band'][i]['frequency'] for i in range(nb)])


xt = []
xl = []
for i in range(nq):
    try:
        xl.append(yam['phonon'][i]['label'])
        xt.append(yam['phonon'][i]['distance'])
    except:
        continue

THz_to_cminv = 33.35641
freqs = transpose(freqs)*THz_to_cminv

figure()
plot(qp, freqs, color = 'k')
xlim(qp[0],qp[-1])
ylim(freqs.min()*1.05, freqs.max()*1.05)
xticks(xt,xl)
ylabel('Frequency (cm$^{-1}$)')
tight_layout()
show()
