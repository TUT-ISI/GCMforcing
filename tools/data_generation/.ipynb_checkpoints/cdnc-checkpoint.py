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

def lnS_crit_insol_over_A(pb,pd):

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
    if zgamma_m >= 0.:

        tmp = (0.5*np.abs(zgamma_m))**inv3

    else:

        tmp = -(0.5*np.abs(zgamma_m))**inv3
        
    zone_thirds = (0.5*zgamma_p)**inv3        + tmp #sign((0.5*np.abs(zgamma_m))**inv3,zgamma_m)
    
    #raising zgamma_m to the second power first takes care of negative values here
    ztwo_thirds = (((0.5*zgamma_p)**(2.))**inv3)+(((0.5*zgamma_m)**(2.))**inv3)
    
    res = (2.*zb_third**2+zb_third*zone_thirds+ztwo_thirds)/\
        (9.*zb_third**3+6.*zb_third**2*zone_thirds+3.*zb_third*ztwo_thirds+pd) 
    
    return res

def cloud_activation(pnaero, pvols, ptm1, papm1, pqm1, pw, nw):#, pwpdf):
    
    zeps=sys.float_info.epsilon
    
    # Raoult effect [m3/mol] NOTE bb must be multiplied by the number of moles of # solute
    zbb = 6.*(molarweight['WAT'])/ \
        (np.pi*density['WAT'])   
    pnact=dict()
    pfracn=dict()
    pcdncact=dict()
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
    
    for b in bins:

        zvols_sol[b]=np.zeros(dimensions)
        zvols_insol[b]=np.zeros(dimensions)
        zvols_tot[b]=np.zeros(dimensions)
        pnact[b]=np.zeros(dimensions)
        pfracn[b]=np.zeros(dimensions)
        pcdncact=np.zeros(dimensions)
        
        
        for spec in specs[2:]:
            # species in a bin 
            key = '{}_{}'.format(spec,b)
            zvols_sol[b]+=pvols[key]*soluble[spec]
            zvols_insol[b]+=pvols[key]*(1-soluble[spec])
            zvols_tot[b]+=zvols_sol[b]+zvols_insol[b]

    for kk in range(dimensions[3]):
    
        for jj in range(dimensions[2]):
            
            for ii in range(dimensions[1]):
                
                for time in range(dimensions[0]):
                    
                    zaa = 4.*(molarweight['WAT'])*surft_w/ \
                        (argas*(density['WAT'])*ptm1[time,ii,jj,kk]) # Kelvin effect [m]

                    if(pw[time,ii,jj,kk] > zeps and pqm1[time,ii,jj,kk] > zeps and ptm1[time,ii,jj,kk] > cthomi):  
                        
                        zntot  = 0.
                        zsum1  = 0.
                        
                        #-- subrange 1a -- neglected
                        
                        #-- subrange 2a + 2b

                        for b in bins[3:]:
                            
                            numkey = '{}_{}'.format(specs[0],b)
                            
                            if pnaero[numkey][time,ii,jj,kk] > nlim:

                                zns[b]=0.
                                #-- number of moles of solute in one particle [mol]
                                for spec in specs[2:]:
                                    # name of the chemical compound
                                    key = '{}_{}'.format(spec,b)
                                    
                                    zns[b]+=nu[spec]*pvols[key][time,ii,jj,kk] * density[spec]/\
                                        molarweight[spec]/pnaero[numkey][time,ii,jj,kk]
                                    
                                #-- critical supersaturation, Koehler equation
                                zcc = np.sqrt(3.*zbb*zns[b]/zaa)
                                zdd = zvols_insol[b][time,ii,jj,kk]*6./(np.pi*pnaero[numkey][time,ii,jj,kk])
                                tmp = lnS_crit_insol_over_A(zcc, zdd)
                                psc[b] = np.exp(zaa*tmp)-1.

                                #-- sums in equation (8), part 3
                                zntot+=pnaero[numkey][time,ii,jj,kk]
                                zsum1+=pnaero[numkey][time,ii,jj,kk]/psc[b]**(2./3.)
                                
                        if(zntot < nlim  or  zsum1 < nlim):
                            continue
            
                        #-- latent heat of evaporation [J/kg]
                        zevap  = 2.501e6-2370.*(ptm1[time,ii,jj,kk]-273.15)
                        
                        #-- saturation vapor pressure of water [Pa]
                        #   Seinfeld \ Pandis (1.10)
                        za1    = 1.-(373.15/ptm1[time,ii,jj,kk])
                        zps    = p0sl_bg*np.exp(13.3185*za1-1.976*za1**2-0.6445*za1**3-0.1299*za1**4)     
                        
                        #-- part 1, eq (11)
                        zalpha = grav*(molarweight['WAT'])*zevap/(cpd*argas*ptm1[time,ii,jj,kk]**2)-\
                            grav*mair/(argas*ptm1[time,ii,jj,kk])
                        
                        #-- part 1, eq (12)
                        zgamma = argas*ptm1[time,ii,jj,kk]/(zps*molarweight['WAT']) \
                            + molarweight['WAT']*zevap**2/(cpd*papm1[time,ii,jj,kk]*mair*ptm1[time,ii,jj,kk])
                        
                        #-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
                        #  Eq (17.61) in second edition
                        zdv= 1.e-4 * (0.211*p0sl_bg/papm1[time,ii,jj,kk]) * ((ptm1[time,ii,jj,kk]/tmelt)**1.94)
                        
                        #-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
                        # Eq (17.71) in second edition
                        zka= 1.e-3 * (4.39 + 0.071 * ptm1[time,ii,jj,kk])
                        
                        #-- growth coefficient, part 1, eq (16)
                        #-- (note: here uncorrected diffusivities and conductivities are used
                        #    based on personal communication with H. Abdul-Razzak, 2007)
                        zgc = 1./(density['WAT']*argas*ptm1[time,ii,jj,kk]/                \
                                  (zps*zdv*molarweight['WAT']) +                      \
                                  zevap*density['WAT']/(zka*ptm1[time,ii,jj,kk]) *         \
                                  (zevap*molarweight['WAT']/(ptm1[time,ii,jj,kk]*argas)-1.))
            
                        #-- effective critical supersaturation: part 3, eq (8)
                        zseff = (zntot/zsum1)**(3./2.)
                        
                        #-- part 3, equation (5)
                        ztheta = ((zalpha*pw[time,ii,jj,kk]/zgc)**(3./2.))/\
                            (2.*np.pi*density['WAT']*zgamma*zntot)
                        
                        #-- part 3, equation (6)
                        zkhi = (2./3.)*zaa*np.sqrt(zalpha*pw[time,ii,jj,kk]/zgc)
                        
                        #-- maximum supersaturation of the air parcel: part 3, equation (9)
                        zsum2 = 0.5*(zkhi/ztheta)**(3./2.)     \
                            + ((zseff**2)/(ztheta+3.*zkhi))**(3./4.)
                        
                        if (zsum2 > zeps) :
                            psmax = zseff / np.sqrt(zsum2)
                        else:
                            psmax=0.
                                
                        bin_index = -1

                        sum_nact = 0.

                        for b in bins[3:]:

                            zfrac[b]=0.
                            
                            numkey = '{}_{}'.format(specs[0],b)
                            
                            bin_index+=1
                            if (pnaero[numkey][time,ii,jj,kk] > nlim  and  zns[b] > 0.):
                                
                                #-- moles of solute in particle at the upper bound of the bin
                                znshi = zns[b]*pnaero[numkey][time,ii,jj,kk]*vhilim[bin_index]/zvols_tot[b][time,ii,jj,kk]
                                
                                #-- lower bound of critical supersaturation
                                zcc = np.sqrt(3.*zbb*znshi/zaa)
                                zdd = (6./np.pi)*vhilim[bin_index]*zvols_insol[b][time,ii,jj,kk]/zvols_tot[b][time,ii,jj,kk]
                                tmp = lnS_crit_insol_over_A(zcc, zdd)
                                sil = np.exp(zaa*tmp)-1.

                                if psmax < sil :
                                    continue
                                
                                #-- moles of solute at the lower bound of the bin:
                                znslo = zns[b]*pnaero[numkey][time,ii,jj,kk]*vlolim[bin_index]/zvols_tot[b][time,ii,jj,kk]
                                
                                #-- upper bound of critical supersaturation
                                zcc = np.sqrt(3.*zbb*znslo/zaa)
                                zdd = (6./np.pi)*vlolim[bin_index]*zvols_insol[b][time,ii,jj,kk]/zvols_tot[b][time,ii,jj,kk]
                                tmp = lnS_crit_insol_over_A(zcc, zdd)
                                siu = np.exp(zaa*tmp)-1.
                                #-- fraction of activated particles in bin, eq (13), part 3
                                zfrac[b] = np.minimum(1.,np.log(psmax/sil)/np.log(siu/sil))
                                pfracn[b][time,ii,jj,kk] = zfrac[b]
                                pnact[b][time,ii,jj,kk] = pfracn[b][time,ii,jj,kk]*pnaero[numkey][time,ii,jj,kk]
                                pcdncact[time,ii,jj,kk]+= pnact[b][time,ii,jj,kk]#*pwpdf[jw]
                                if ii==47:
                                    print(b,pfracn[b][time,ii,jj,kk])
                        
    return pcdncact
                
# end def cloud_activation

  
