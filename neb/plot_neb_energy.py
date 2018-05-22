from pylab import *
rcParams.update({'font.size': 48, 'text.usetex': True})


nebstep = arange(0,6)
energy  = array([-80.51356837 ,
				 -78.67600771 ,
				 -80.50738854 ,
				 -80.50740238 ,
				 -80.50736986 ,
				 -80.51266666 ])
print("Barrier height = %f"%(max(energy) - min(energy)))
figure()
plot(nebstep, energy)
xlabel("NEB step")
ylabel("Total Energy (eV)")
xticks(nebstep)
tight_layout()
show()
