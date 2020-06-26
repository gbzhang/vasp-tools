from pylab import *
from scipy.optimize import leastsq
import scipy.constants as constants
import mendeleev
import os, sys
rcParams.update({'font.size':48, 'text.usetex': True})


def BirchMurnaghan(parameters, vol, order):
    """ 
    Summary:  
      Fitting routine for Birch Murnaghan of orders 3 to 10
    Inputs:
      parameters: array of coefficients: [V0, E0, B0, BP, c4, c5, ... ]
      vol:        array of volume points
      order:      order of BM fit 
    Outputs:
      E: energy     resulting energy array
      P: pressure   resulting pressure array
    """
    V0 = parameters[0]
    E0 = parameters[1]
    B0 = parameters[2]
    BP = parameters[3]

    f = ( (V0/vol)**(2./3.) - 1 )/2

    if order == 3:
        E = E0 + V0*(9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 )
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 )
        
    elif order == 4:
        c4 = parameters[4];
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 )
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 )
        
    elif order == 5:
        c4=parameters[4]; c5=parameters[5];
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 + B0*c5/120*f**5 )
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 + B0*c5/24*f**4 )

    elif order == 6:
        c4=parameters[4]; c5=parameters[5]; c6=parameters[6];
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 + B0*c5/120*f**5 + B0*c6/720*f**6)
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 + B0*c5/24*f**4 + B0*c6/120*f**5 )

    elif order == 7:
        c4=parameters[4]; c5=parameters[5]; c6=parameters[6]; c7=parameters[7]; 
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 + B0*c5/120*f**5 + B0*c6/720*f**6 + B0*c7/5040*f**7)
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 + B0*c5/24*f**4 + B0*c6/120*f**5 + B0*c7/720*f**6 )

    elif order == 8:
        c4=parameters[4]; c5=parameters[5]; c6=parameters[6]; c7=parameters[7]; c8=parameters[8];
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 + B0*c5/120*f**5 + B0*c6/720*f**6 + B0*c7/5040*f**7 + B0*c8/40320*f**8)
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 + B0*c5/24*f**4 + B0*c6/120*f**5 + B0*c7/720*f**6 +
                                     B0*c8/5040*f**7)

    elif order == 9:
        c4=parameters[4]; c5=parameters[5]; c6=parameters[6]; c7=parameters[7]; c8=parameters[8]; c9=parameters[9];
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 + B0*c5/120*f**5 + B0*c6/120*f**6 +
                      B0*c7/5040*f**7 + B0*c8/40320*f**8 + B0*c9/362880*f**9 )
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 + B0*c5/24*f**4 + B0*c6/120*f**5 + B0*c7/720*f**6 +
                                     B0*c8/5040*f**7 + B0*c9/40320*f**8 )

    elif order == 10:
        c4=parameters[4]; c5=parameters[5]; c6=parameters[6]; c7=parameters[7]; c8=parameters[8]; c9=parameters[9]; c10=parameters[10];
        E = E0 + V0*( 9.*B0/2*f**2 + 27.*B0*(BP-4.)/6*f**3 + B0*c4/24*f**4 + B0*c5/120*f**5 + B0*c6/120*f**6 +
                      B0*c7/5040*f**7 + B0*c8/40320*f**8 + B0*c9/362880*f**9 + B0*c10/3628800*f**10 )
        P = (V0/vol)**(5./3.)/3. * ( 9.*B0*f + 27.*B0*(BP-4.)/2*f**2 + B0*c4/6*f**3 + B0*c5/24*f**4 + B0*c6/120*f**5 + B0*c7/720*f**6 +
                                     B0*c8/5040*f**7 + B0*c9/40320*f**8 + B0*c10/362880*f**9)

    else:
        print('Error: only orders 3-10 are allowed')
        exit(1)
        
    return E, P


def initialguess(e, v, vfit, order):
    """
    Summary:
      provides an initial guess for the parameters based on a simple quadratic: 
         y = ax^2 + bx + c
    Inputs:
      e:      array of calculated energy points
      v:      array of calculated volume points
      vfit:   dense array of volume points to interpolate between points in v
      order:  order of the BM fitting routine to be used
    Outputs:
      x0:     array of initial guess parameters for a given order, [V0, E0, B0, BP, c4, ... ]
    """
    a,b,c = polyfit(v,e,2)
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4.
    if order == 3:
        x0 = [v0, e0, b0, bP]
    elif order == 4:
        x0 = [v0, e0, b0, bP, 0.1]
    elif order == 5:
        x0 = [v0, e0, b0, bP, 0.1, 0.1]
    elif order == 6:
        x0 = [v0, e0, b0, bP, 0.1, 0.1, 0.1]
    elif order == 7:
        x0 = [v0, e0, b0, bP, 0.1, 0.1, 0.1, 0.1]
    elif order == 8:
        x0 = [v0, e0, b0, bP, 0.1, 0.1, 0.1, 0.1, 0.1]
    elif order == 9:
        x0 = [v0, e0, b0, bP, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    elif order == 10:
        x0 = [v0, e0, b0, bP, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    return x0


def objective(pars,y,x,order):
    """
    Summary: objective function for passing to the scipy leastsq routine
    """
    E, P = BirchMurnaghan(pars, x, order)
    err =  y - E

    return err


def getfit(v, e, vfit, order):
    """
    Summary: 
      Handling routine to compute E and P from given v, e input and given order
    Inputs:
      v:     array of calculated volume points
      e:     array of calculated energy points
      vfit:  array of interpolating volume points to do the fitting
      order: BM fit order
    Outputs:
      E:        calculated energy array at vfit points
      P:        calculated pressure array at vfit points
      murnpars: minimization parametrs from BM fitting
    """
    pconv = constants.e*1e21 # convert pressure from eV/Ang^3 to GPa
    x0 = initialguess(e,v,vfit,order)
    murnpars, ier = leastsq(objective, x0, args=(e,v,order),ftol=1e-12,xtol=1e-12)
    E, P = BirchMurnaghan(murnpars, vfit, order)
    murnpars[2] *= pconv
    print('%16s'%order, *["%16.5f"%i for i in murnpars[:4]])

    return E, P, murnpars

if __name__ == '__main__':
    w     = mendeleev.Sn.mass  # atomic mass g/mol
    eV    = constants.e        # electron charge/ eV to J
    NA    = constants.N_A      # Avogadro's number
    pconv = eV*1e21            # convert pressure from eV/Ang^3 to GPa

    data = genfromtxt('energies.dat')
    v = data[:,0]
    e = data[:,1]

    # capability for printing density in g/cm^3
    # print('#%15s %16s' %('volume [cc/g]','energy (kJ/g)'))
    # [print("%16.8f %16.8f"%(v[i], e[i])) for i in range(len(v))]
    print('%16s %16s %16s %16s %16s'%('order', 'V0 (Ang^3/atom)', 'E0 (eV/atom)', 'B0 (GPa)', 'BP'))
    
    V = linspace(v[0],v[-1],1000)
    E3,   P3,  mp3 = getfit(v,e,V,3)
    E4,   P4,  mp4 = getfit(v,e,V,4)
    E5,   P5,  mp5 = getfit(v,e,V,5)
    E6,   P6,  mp6 = getfit(v,e,V,6)
    E7,   P7,  mp7 = getfit(v,e,V,7)
    E8,   P8,  mp8 = getfit(v,e,V,8)
    E9,   P9,  mp9 = getfit(v,e,V,9)
    E10, P10, mp10 = getfit(v,e,V,10)

    # Plot energy
    figure()
    plot(v,  e, 'o', color='C0')
    plot(V, E3, '-', label='order 3')
    plot(V, E4, '-', label='order 4')
    plot(V, E5, '-', label='order 5')
    plot(V, E6, '-', label='order 6')
    plot(V, E7, '-', label='order 7')
    plot(V, E8, '-', label='order 8')
    plot(V, E9, '-', label='order 9')
    plot(V, E10, '-', label='order 10')
    legend(fontsize=24)
    xlabel('Volume ($\mathrm{\AA}^3$/atom)')
    ylabel('Energy (eV/atom)')
    #xlim(20,33)
    #ylim(-4,-3.5)
    tight_layout()
    show()

    # Plot pressure
    figure()
    plot(V, P3*pconv, '-', label='order 3')
    plot(V, P4*pconv, '-', label='order 4')
    plot(V, P5*pconv, '-', label='order 5')
    plot(V, P6*pconv, '-', label='order 6')
    plot(V, P7*pconv, '-', label='order 7')
    plot(V, P8*pconv, '-', label='order 8')
    plot(V, P9*pconv, '-', label='order 9')
    plot(V, P10*pconv, '-', label='order 10')
    legend(fontsize=24)
    xlabel('Volume ($\mathrm{\AA}^3$/atom)')
    ylabel('Pressure (GPa)')
    #xlim(20,33)
    #ylim(-4,-3.5)
    tight_layout()
    show()
