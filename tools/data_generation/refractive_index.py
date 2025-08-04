from salsa_parameters import bins, specs
from parameters import refrac
import numpy as np
from ref_ind_parameters import import_refrac

def refractive_index(*args):
    """
    by  Atte Laakso / UEF

    Reads data for refractive indices based on given arguments
    """
    
    # if only one parameter is given, get refractive indices from SALSA parameters
    if len(args)==1:

        # particle volumes
        pvols=args[0]
        # temporary variables
        nrpart=dict()
        nipart=dict()
        # real part of the refractive index 
        nr=dict()
        # imaginary part of the refractive index 
        ni=dict()

        # dimensions of the dataset
        dimensions=pvols['SO4_1a1'].shape

        # volume data
        vols_tot=np.zeros(dimensions)
        vtot=vols_tot
        
        # iterate over different bins
        for b in bins:
            # initialization of variables
            nrpart[b]=np.zeros(dimensions)
            nipart[b]=np.zeros(dimensions)
            vols_tot=np.zeros(dimensions)
            for spec in specs[1:]:
                # name of the chemical compound 
                key = '{}_{}'.format(spec,b)
                # total particle volume
                vols_tot+=np.absolute(pvols[key])
                # volume weighting of refractive index
                nrpart[b]+=np.absolute(pvols[key])*refrac[spec].real
                nipart[b]+=np.absolute(pvols[key])*refrac[spec].imag
        
            # final refractive index
            # (set limits for empty grid cells)
            vtot=np.maximum(1e-30,vols_tot)
            # 
            nr[b]=np.maximum(nrpart[b]/vtot,1.33)
            # 
            ni[b]=np.maximum(nipart[b]/vtot,1.e-9)

    
    # if two parameters are given, get refractive indices from model's data
    elif len(args) == 2:
        print("Reading in refractive indices from the model data.")
        # Get model name
        model_name=args[1]
        # get pvols
        pvols=args[0]
        
        # temporary variables
        nrpart=dict()
        nipart=dict()
        # real part of the refractive index 
        nr=dict()
        # imaginary part of the refractive index 
        ni=dict()
        
        dimensions=pvols['SO4_1a1'].shape
        
        vols_tot=np.zeros(dimensions)
        vtot=vols_tot

        #get refractive indices of the model
        mrefrac=import_refrac(model_name)
        
        for b in bins:
            # initialization of variables
            nrpart[b]=np.zeros(dimensions)
            nipart[b]=np.zeros(dimensions)
            vols_tot=np.zeros(dimensions)
            
            for spec in specs[1:]:
                # The refractive indices data doesn't have all the parameters given as specific as in salsa_parameters
                #
                #
                #
                # Alter this part if needed!!!
                # 
                #
                #
                # Read in the best possible value for these substances
                if spec in ("VBS1" , "VBS10" , "IEPOX" , "Glyx" , "OC") :
                    
                    # searching name for the best corresponding compound (organic aerosols)
                    a_name='OA'                   
                    # name of the chemical compound 
                    key = '{}_{}'.format(spec,b)
                    # total particle volume
                    vols_tot+=np.absolute(pvols[key])
                    
                    # volume weighting of refractive index
                    nrpart[b]+=np.absolute(pvols[key])*mrefrac[a_name].real
                    nipart[b]+=np.absolute(pvols[key])*mrefrac[a_name].imag
                    
                else:
                    # name of the chemical compound 
                    key = '{}_{}'.format(spec,b)
                    
                    # total particle volume
                    vols_tot+=np.absolute(pvols[key])
                    
                    # volume weighting of refractive index
                    nrpart[b]+=np.absolute(pvols[key])*mrefrac[spec].real
                    nipart[b]+=np.absolute(pvols[key])*mrefrac[spec].imag
            
            # final refractive index
            # (set limits for empty grid cells)
            vtot=np.maximum(1e-30,vols_tot)
            # 
            nr[b]=np.maximum(nrpart[b]/vtot,1.33)
            # 
            ni[b]=np.maximum(nipart[b]/vtot,1.e-9)
        
    else:
        print("Reafracive_index accepts only one or two parameters.")

    return nr, ni
