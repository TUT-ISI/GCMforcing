# fraction of organic compounds hygroscopicity compared to an ideal solution
epsoc=0.15
# molar weight [kg mol-1]
molarweight={
    "WAT"   : 18.016e-3,
    "SO4"   : 98.08e-3,
    "OC"    : 180.e-3,
    "SS"    : 58.44e-3,
    "BC"    : 12.e-3,
    "DU"    : 100.e-3,
    "VBS1"  : 186.e-3,
    "VBS10" : 186.e-3,
    "IEPOX" : 118.e-3,
    "Glyx"  : 58.e-3,
    }

# binary molality of soluble species as a function of water activity
zbinmol={
    "SO4" : lambda zaw: 1.1065495e2-3.6759197e2*zaw+5.0462934e2*zaw**2-3.1543839e2*zaw**3+6.770824e1*zaw**4,
    "OC"  : lambda zaw: (1.0/(zaw*molarweight["OC"])-1.0/molarweight["OC"])/epsoc,
    "SS"  : lambda zaw: 5.875248e1-1.8781997e2*zaw+2.7211377e2*zaw**2-1.8458287e2*zaw**3+4.153689e1*zaw**4,
     }

# density [kg m-3]
density={
    "WAT"   : 1000.,
    "SO4"   : 1830.,
    "OC"    : 2000.,
    "SS"    : 2165.,
    "BC"    : 2000.,
    "DU"    : 2650.,
    "VBS1"  : 1320.,
    "VBS10" : 1320.,
    "IEPOX" : 1320.,
    "Glyx"  : 1320.,
    }

# is the chemical compound soluble (1) or not (2)
soluble={
    "WAT"   : 0.,
    "SO4"   : 1.,
    "OC"    : 0.,
    "SS"    : 1.,
    "BC"    : 0.,
    "DU"    : 0.,
    "VBS1"  : 1.,
    "VBS10" : 1.,
    "IEPOX" : 1.,
    "Glyx"  : 1.,
    }

# dissociation constant of the chemical compound
nu={
    "WAT"   : 0.,
    "SO4"   : 3.,
    "OC"    : 1.,
    "SS"    : 2.,
    "BC"    : 0.,
    "DU"    : 0.,
    "VBS1"  : 1.,
    "VBS10" : 1.,
    "IEPOX" : 1.,
    "Glyx"  : 1.,
    }

# surface tension [J m-2]
surft_w = 0.073 # @ 293 K

# Gas constant
argas = 8.314 # J K-1 mol

# Avogadro's constant
avog=6.0221e+3

# Boltzmann constant [J K-1]
boltzmann = 1.3807e-23

# Molar weight of air [mol kg-1]
mair = 28.97e-3

# Standard pressure [Pa]
p0sl_bg         = 101325.

# Melting temperature of ice/snow
tmelt = 273.15

# Minimum temperature for mixed clouds
cthomi  = tmelt-35.0

# Specific heat at constant pressure
cpd   = 1004.64

# Gravitational acceleration [m s-2]
grav  = 9.80665

# Low limit of number concentration
nlim = 1.0

# refractive indices of chemical compounds
refrac={
    "WAT"   : complex(1.335, 2.8e-9),
    "SO4"   : complex(1.432, 1.e-9),
    "OC"    : complex(1.53, 5.5e-3),
    "SS"    : complex(1.45, 1.e-8),
    "BC"    : complex(1.85, 7.1e-1),
    "DU"    : complex(1.45, 1.e-3),
    "VBS1"  : complex(1.53, 5.5e-3),
    "VBS10" : complex(1.53, 5.5e-3),
    "IEPOX" : complex(1.53, 5.5e-3),
    "Glyx"  : complex(1.53, 5.5e-3),
    }


