import numpy as np

def asymmetry(pxx,pnr,pni,lookup):

    x0_min=np.array([0.001, 0.16, 5.0e-6, 0.0015])     #min Mie parameter
    x0_max=np.array([25.0, 210.0, 3.0, 17.0])          #max Mie parameter
    nr_min=np.array([1.33, 1.33, 1.0, 1.0])            #min ref. ind. real
    nr_max=np.array([2.00, 2.00, 3.0, 3.0])            #max ref. ind. real
    ni_min=np.array([1.0e-9,  1.0e-9, 1.0e-9, 1.0e-9]) #min ref. ind. imag.
    ni_max=np.array([1.0, 1.0, 2.0, 2.0])              #max ref. ind. imag.

    log_x0_min=np.log(x0_min)
    log_x0_max=np.log(x0_max)
    log_nr_min=np.log(nr_min)
    log_nr_max=np.log(nr_max)
    log_ni_min=np.log(ni_min)
    log_ni_max=np.log(ni_max)
    
    # the dimensions of the look-up tables (totally unnecessary):
    Ndismax =np.array([100, 100, 100, 100]) 
    Nnrmax  =np.array([100, 100, 100, 100])
    Nnimax  =np.array([200, 200, 200, 200])

    # linearized value increment -- for some reason not implemented
    # for size parameter:
    inc_nr=(nr_max-nr_min)/Nnrmax
    inc_ni=(log_ni_max-log_ni_min)/Nnimax

    # check, which table to use depending on the size parameter
    ktable = np.where(pxx >= x0_max[0], 1, 0)

    # indices in the look-up-table that correspond to size parameters pxx
    Ndis=np.array((np.log(pxx)-log_x0_min[0]) /
                  (log_x0_max[0]-log_x0_min[0])*Ndismax[0]*(1-ktable)
                  +
                  (np.log(pxx)-log_x0_min[1]) /
                  (log_x0_max[1]-log_x0_min[1])*Ndismax[1]*(ktable), np.int32)
          
    Ndis=np.minimum(Ndismax[0]-1,np.maximum(0,Ndis))
          
    # indices in the look-up-table that correspo nd to
    # the real part of refractive index pnr
    Nnr=np.array((pnr-nr_min[0])/inc_nr[0]*(1-ktable)
                 +
                 (pnr-nr_min[1])/inc_nr[1]*(ktable), np.int32)
          
    Nnr=np.minimum(Nnrmax[0]-1,np.maximum(0,Nnr))

    # indices in the look-up-table that correspond to
    # the imaginary part of refractive index pni
    Nni=np.array((np.log(pni)-log_ni_min[0])/inc_ni[0]*(1-ktable)
                 +
                 (np.log(pni)-log_ni_min[1])/inc_ni[1]*(ktable), np.int32)

    Nni=np.minimum(Nnimax[0]-1,np.maximum(0,Nni))

    # asym_n = asymmetry factor in table n
    asym_1='asym_1'
    asym_2='asym_2'

    # read 3D asym_n from the lut
    lut1=(np.array(lookup.variables['asym_1']))
    lut2=(np.array(lookup.variables['asym_2']))

    # retrieve the ASY values corresponding to pxx, pnr, and pni 
    ASY=lut1[Ndis,Nni,Nnr]*(1-ktable)+lut2[Ndis,Nni,Nnr]*ktable 

    return ASY

# do not run main if imported
if __name__ == "__main__":
    main()
