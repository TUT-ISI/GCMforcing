
##################################################################################
## \brief
## Implementation of Abdul-Razzak & Ghan activation scheme for SALSA.
##
## \author Harri Kokkola (FMI)
##
## \responsible_coder
## Thomas Kuehn -- thomas.h.kuhn@uef.fi
## Harri Kokkola -- harri.kokkola@fmi.fi
##
## \revision_history
## T. Anttila (FMI)     2007
## H. Kokkola (FMI)     2007
## A.-I. Partanen (FMI) 2007
## T. Kuehn (UEF)       2015, 2018
## H. Kokkola (FMI)     2019
##
## \limitations
## None
##
## \details
## Purpose: Calculates the number of activated cloud 
## droplets according to parameterizations by:
##
## Abdul-Razzak et al: "A parameterization of aerosol activation - 
##                      3. Sectional representation"
##                      J. Geophys. Res. 107, 10.1029/2001JD000483, 2002. 
##                      [Part 3]
##
## Abdul Razzak et al: "A parameterization of aerosol activation - 
##                      1. Single aerosol type"
##                      J. Geophys. Res. 103, 6123-6130, 1998. 
##                      [Part 1]
##
## \copyright
## The script is licensed under the Apache License 2.0
##
##################################################################################
import numba as nb
import sys
from parameters import *
from salsa_parameters import *
import numpy as np

def lnS_crit_insol_over_A_dp(pb,pd):

    # -----------------------------------------------------------------------
    #
    # def lnS_crit_insol_over_A
    #
    # Computes ln(S_crit) for the case that the aerosol particle containes an 
    # insoluble core. The equations where derived starting from Eq. (17.38) in
    # Seinfeld \ Pandis, 2nd edition, essentially solving a cubic equation of 
    # the form 
    #
    # x**3 + bx**2 + d = 0
    #
    # where b = SQRT(3B/A) and d = du**3. Here A and B denote the Kelvin and
    # Raoult effects, respectively, and du is the diameter of the insoluble 
    # core. The result is normalised with the parameter A.
    #
    # Seinfeld \ Pandis didn't provide a solution to the problem and I
    # didn't find a reference either, so I did the calculations myself. If
    # anybody knows of a reference, please let me know. In case of any doubt,
    # I can provide the calculations and a python code that tests the results.
    #  
    # Authors:
    # --------
    # Thomas Kuehn (thk), UEF    6/2018 --
    #
    # lnS_crit_insol_over_A is called from cloud_activation in mo_ham_salsa_cloud
    #
    # -----------------------------------------------------------------------
    
    inv3 = 1./3.
    
    zb_third = inv3*pb
    
    zgamma_p = 2.*zb_third**3+pd+np.sqrt(4.*zb_third**3*pd+pd**2)
    zgamma_m = 2.*zb_third**3+pd-np.sqrt(4.*zb_third**3*pd+pd**2)
    
    #-->eehol: when number of moles of solute is low the zgamma_m becomes negative which messes up the cubic root
    #taking the cubic root from absolute value and putting the original sign of zgamma_m takes care of this
    zone_thirds = (0.5*zgamma_p)**inv3 + np.sign(zgamma_m) * (0.5*np.abs(zgamma_m))**inv3
    
    #raising zgamma_m to the second power first takes care of negative values here
    ztwo_thirds = (((0.5*zgamma_p)**(2.))**inv3)+(((0.5*zgamma_m)**(2.))**inv3)
    
    res = (2.*zb_third**2+zb_third*zone_thirds+ztwo_thirds)/\
        (9.*zb_third**3+6.*zb_third**2*zone_thirds+3.*zb_third*ztwo_thirds+pd) 
    
    return res

def cloud_activation_dp(pnaero, pvols, ptm1, papm1, pqm1, pw, nw, psmax, psmax_fixed, ratio):
    
    zeps=sys.float_info.epsilon
    
    # Raoult effect [m3/mol] NOTE bb must be multiplied by the number of moles of # solute
    zbb = 6.*(molarweight['WAT'])/ \
        (np.pi*density['WAT'])   
    # Kelvin effect [m]
    zaa = 4.*(molarweight['WAT'])*surft_w/ \
        (argas*(density['WAT'])*ptm1)

    pnact=dict()          # number of activated droplets in a bin #/m3
    pfracn=dict()         # fractions of activated droplets in a bin #/3
    psc=dict()
    zfrac=dict()
    prc  = -1.  # this is not calculated here
    
    # Keeping track of what is soluble and what is insoluble
    zvols_sol=dict()
    zvols_insol=dict()
    zvols_tot=dict()

    # initialize
    dimensions=ptm1.shape
    zns=dict()
    znshi=dict()
    znslo=dict()
    zntot=dict()
    pcdncact=dict()
    pcdncact_final=dict()
    pcdncact=np.zeros(dimensions)
    pcdncact_final=pcdncact
    zntot=np.zeros(dimensions)
    zsum1=np.zeros(dimensions)
    znshi=np.zeros(dimensions)
    znslo=np.zeros(dimensions)

    if psmax_fixed==0:
        psmax=np.ones(dimensions)

    zseff=np.ones(dimensions)

    #cloud_mask=np.where((pw>zeps)&(pqm1>zeps)&(ptm1>cthomi))

    for b in bins:

        zvols_sol[b]=np.zeros(dimensions)
        zvols_insol[b]=np.zeros(dimensions)
        zvols_tot[b]=np.zeros(dimensions)
        pnact[b]=np.zeros(dimensions)
        pfracn[b]=np.zeros(dimensions)
        zns[b]=np.zeros(dimensions)
        
        for spec in specs[2:]:
            # species in a bin 
            key = '{}_{}'.format(spec,b)
            numkey = '{}_{}'.format(specs[0],b)
            # sum of volumes of all soluble compounds in bin b
            zvols_sol[b]+=pvols[key]*soluble[spec]
            # sum of volumes of all insoluble compounds
            zvols_insol[b]+=pvols[key]*(1.-soluble[spec])
            # total volume of all compounds
            zvols_tot[b]+=pvols[key]*soluble[spec]+pvols[key]*(1.-soluble[spec])
            # number of moles of soluble compounds in a particles
            zns[b]+=nu[spec]*pvols[key] * density[spec]/\
                molarweight[spec]/np.maximum(pnaero[numkey],nlim)
                        
        if psmax_fixed==0:
            
            #-- critical supersaturation, Koehler equation
            zcc = np.sqrt(3.*zbb*zns[b]/zaa)
            zdd = zvols_insol[b]*6./(np.pi*np.maximum(pnaero[numkey],nlim))
            tmp = lnS_crit_insol_over_A_dp(zcc, zdd)
            psc[b] = np.exp(zaa*tmp)-1.
            #psc[b] = np.exp(4.*zaa**3/(27.*zbb)/zns[b]) - 1.
        
            if b[0]!='1': # subregion 1 is neglected
            
                #-- sums in equation (8), part 3
                zntot+=pnaero[numkey] 
                zsum1+=pnaero[numkey]/psc[b]**(2./3.)

    if psmax_fixed==0:
        
        zevap  = 2.501e6-2370.*(ptm1-273.15)
                        
        #-- saturation vapor pressure of water [Pa]
        #   Seinfeld \ Pandis (1.10)
        za1    = 1.-(373.15/ptm1)
        zps    = p0sl_bg*np.exp(13.3185*za1-1.976*za1**2-0.6445*za1**3-0.1299*za1**4)     
        
        #-- part 1, eq (11)
        zalpha = grav*(molarweight['WAT'])*zevap/(cpd*argas*ptm1**2)-\
            grav*mair/(argas*ptm1)
                        
        #-- part 1, eq (12)
        zgamma = argas*ptm1/(zps*molarweight['WAT']) \
            + molarweight['WAT']*zevap**2/(cpd*papm1*mair*ptm1)
                        
        #-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
        #  Eq (17.61) in second edition
        zdv= 1.e-4 * (0.211*p0sl_bg/papm1) * ((ptm1/tmelt)**1.94)
                        
        #-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
        # Eq (17.71) in second edition
        zka= 1.e-3 * (4.39 + 0.071 * ptm1)
                        
        #-- growth coefficient, part 1, eq (16)
        #-- (note: here uncorrected diffusivities and conductivities are used
        #    based on personal communication with H. Abdul-Razzak, 2007)
        zgc = 1./(density['WAT']*argas*ptm1/                \
                  (zps*zdv*molarweight['WAT']) +            \
                  zevap*density['WAT']/(zka*ptm1) *         \
                  (zevap*molarweight['WAT']/(ptm1*argas)-1.))
        
        #-- effective critical supersaturation: part 3, eq (8)
        number_mask=np.where(np.logical_and(zntot>nlim,zsum1>nlim))
        zseff[number_mask] = (zntot[number_mask]/zsum1[number_mask])**(3./2.)
        
        #-- part 3, equation (5)
        ztheta = ((zalpha*pw/zgc)**(3./2.))/\
            (2.*np.pi*density['WAT']*zgamma*zntot)
        
        #-- part 3, equation (6)
        zkhi = (2./3.)*zaa*np.sqrt(zalpha*pw/zgc)
                        
        #-- maximum supersaturation of the air parcel: part 3, equation (9)
        zsum2 = 0.5*(zkhi/ztheta)**(3./2.)     \
            + ((zseff**2)/(ztheta+3.*zkhi))**(3./4.)

        #psmax[cloud_mask] = zseff[cloud_mask]/np.sqrt(zsum2[cloud_mask])
        psmax = zseff/np.sqrt(zsum2)
    
    bin_index = 0

    for b in bins:

        zfrac[b]=np.zeros(dimensions)
                            
        numkey = '{}_{}'.format(specs[0],b)

        if b[0]!='1':            
            #-- moles of solute in particle at the upper bound of the bin
            znshi = zns[b]*pnaero[numkey]*vhilim[bin_index]*ratio/zvols_tot[b]
            
            #-- lower bound of critical supersaturation
            zcc = np.sqrt(3.*zbb*znshi/zaa)
            zdd = (6./np.pi)*vhilim[bin_index]*ratio*zvols_insol[b]/zvols_tot[b]
            tmp = lnS_crit_insol_over_A_dp(zcc, zdd)
            sil = np.exp(zaa*tmp)-1.

            #-- moles of solute at the lower bound of the bin:
            znslo = zns[b]*pnaero[numkey]*vlolim[bin_index]*ratio/zvols_tot[b]
            
            #-- upper bound of critical supersaturation
            zcc = np.sqrt(3.*zbb*znslo/zaa)
            zdd = (6./np.pi)*vlolim[bin_index]*ratio*zvols_insol[b]/zvols_tot[b]
            tmp = lnS_crit_insol_over_A_dp(zcc, zdd)
            siu = np.exp(zaa*tmp)-1.

            p_a=np.where((siu>psmax)& (psmax>sil) & (zns[b]>0.0) & (pnaero[numkey]>nlim))

            tmp=np.minimum(np.log(psmax/sil)/np.log(siu/sil),1.0)
            
            zfrac[b][p_a]=np.maximum(tmp[p_a],0.0)
            # particles with upper bound of critical supersaturation lower than psmax all activate
            zfrac[b][np.where(siu<psmax)]=1.
            
        pfracn[b] = zfrac[b]
        pnact[b] = pfracn[b]*pnaero[numkey]
        pcdncact+= pnact[b]#*pwpdf[jw]
                        
        bin_index+=1

#    pcdncact_final=pcdncact[cloud_mask]
    return pcdncact, pnact
                
# end def cloud_activation

  
