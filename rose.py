from scipy.optimize import least_squares, leastsq
from pylab import *
import os, sys
rcParams.update({'font.size':48, 'text.usetex': True})

def rose(parameters, x):
    
    # rV0=0.7; rE0=1.8; rB0=15; rdBdV=3.5;
    V0 = parameters[3]
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    
    rV0 = V0
    rE0 = E0
    rB0 = B0
    rdBdV = BP
    # E = rE0 + 1000.0*(4.0*rB0*rV0/(rdBdV-1.0)**2.0 -
    #                   2.0*rB0*rV0/(rdBdV-1.0)**2.0 * (5.0+3.0*rdBdV*((x/rV0)**(1./3.)-1.0)-
    #                                                   3.0*(x/rV0)**(1./3.))*exp(3./2.*(rdBdV-1.0)*(1.0-(x/rV0)**(1./3.))))

    E = rE0 + (4.0*rB0*rV0/(rdBdV-1.0)**2.0 -
               2.0*rB0*rV0/(rdBdV-1.0)**2.0 * (5.0+3.0*rdBdV*((x/rV0)**(1./3.)-1.0)-
                                               3.0*(x/rV0)**(1./3.))*exp(3./2.*(rdBdV-1.0)*(1.0-(x/rV0)**(1./3.))))

    return E
    
def BirchMurnaghan(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    equation From PRB 28,5480 (1983)
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    
    E = E0 + (9*V0*B0/16)*( BP*((V0/vol)**(2./3) - 1)**3 + ((V0/vol)**(2./3)-1)**2 * (6-4*(V0/vol)**(2./3)) )
    P = (3*B0/2)*( (V0/vol)**(7./3) - (V0/vol)**(5./3) ) * ( 1 + 0.75*(BP-4)*( (V0/vol)**(2./3) - 1) )

    return E, P


def bm3(parameters, x):
    a0 = parameters[0]
    
    
    E = a0 + 1000.0*aV0*(9.0*a2/2 *(((aV0/x)**(2.0/3.0)-1)/2.0)**2 + 27.0*a2*(a3-4.0)/6 *(((aV0/x)**(2.0/3.0)-1)/2.0)**3 )
    return E

def bm5(parameters, vol):
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    
    cV0 = V0
    c0  = E0
    c2  = B0
    c3  = BP
    c4  = parameters[4]
    c5  = parameters[5]
     
    # fs5 = c0 + 1000.0*cV0*(9.0*c2/2 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**2 +
    #                        27.0*c2*(c3-4.0)/6 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**3 +
    #                        c2*c4/24 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**4+
    #                        c2*c5/120 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**5 )

    fs5 = c0 + cV0*(9.0*c2/2 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**2 +
                    27.0*c2*(c3-4.0)/6 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**3 +
                    c2*c4/24 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**4+
                    c2*c5/120 *(((cV0/vol)**(2.0/3.0)-1)/2.0)**5 )

    
    return fs5

def initialguess(e, v, vfit):
    a,b,c = polyfit(v,e,2)
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4.0
    
    x0 = [e0, b0, bP, v0, 0.1, 0.1]
    #x0 = [e0, b0, bP, v0]
    return x0


def objective(pars,y,x):

    #E, P = BirchMurnaghan(pars, x)
    E = bm5(pars, x)
    #E = rose(pars, x)
    err =  y - E
    return err



if __name__ == '__main__':
    conv    = 1.602e-19*1e21 # eV/Ang^3 to GPa
    data = genfromtxt('energies.dat')

    NA = 6.022e23
    eV = 1.602e-19
    w = 118.71

    v = data[:,0]#*NA/w/1e24
    e = data[:,1]#*NA*eV/w/1000

    print('#%15s %16s' %('volume [cc/g]','energy (kJ/g)'))
    [print("%16.8f %16.8f"%(v[i], e[i])) for i in range(len(v))]
    
    vfit = linspace(v[0],v[-1],1000)
    x0 = initialguess(e,v,vfit)
    murnpars, ier = leastsq(objective, x0, args=(e, v))
    E = bm5(murnpars, vfit)
    print('V0 = ', murnpars[3], ' cm^3/g')
    print('p0 = ', 1/murnpars[3], ' g/cm^3')
    print('B0 = ', murnpars[1]*conv, ' GPa')
    print('E0 = ', murnpars[0], ' kJ/g')

    plot(v,   e, 'o', color='C0')
    plot(vfit, E, '-', color='C0')
    xlabel('Volume ($\mathrm{\AA}^3$)')
    ylabel('Energy (eV)')
    #xlim(15,33)
    #ylim(-5,5)
    tight_layout()
    show()
