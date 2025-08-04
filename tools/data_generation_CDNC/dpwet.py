from parameters import *
from salsa_parameters import *
import numpy as np

def wet_diameter(pnaero,pvols,aw,model,relhum_fixed):
    """
    by  Harri Kokkola / FMI

    Calculates the wet particle size.
    
    Adapted from 'mo_ham_salsa_properties.f90'
    
    """    
    
    dimensions=aw.shape

    # set maximum value of the water activity to a given value

    zaw=aw

    #    if(model == 'CAM-OSLO'):
    zaw=np.minimum(0.95,np.maximum(0.37,aw))
#    zaw=0.37
    # wet diameter [m]
    dwet = dict()
    # dry diameter [m]
    ddry = dict()

    if relhum_fixed > 0:

        zaw=relhum_fixed

    for b in bins:

        # initialize variables
        vols_tot=0.0
        hygr_tot=0.0
        for spec in specs[2:]:
            # name of the chemical compound
            key = '{}_{}'.format(spec,b)
            # total particle volume
            vols_tot+=pvols[key]
            # calculate water taken up by soluble compounds
            if spec in zbinmol.keys():
                hygr_tot+=pvols[key]*density[spec]/molarweight[spec]/zbinmol[spec](zaw)
        # name for the number concentration in a bin (e.g. NUM_1a2)
        numkey = '{}_{}'.format(specs[0],b)
        # dry diameter [m]
        ddry[b]=(6.0/np.pi*vols_tot/pnaero[numkey])**(1.0/3.0)
        # initialize wet diameter 
        dwet[b]=ddry[b]
        # calculate wet diameter

        # no hydration if RH=0%
        if relhum_fixed == 0:

            pvols['WAT']=0.

        # for other cases, calculate water according to RH
        else:
            
            dwet[b]=(hygr_tot/pnaero[numkey]/density['WAT']*6/np.pi+vols_tot/pnaero[numkey]*6.0/np.pi)**(1.0/3.0)

            pvols['WAT']=hygr_tot/density['WAT']

    return dwet, ddry
