"""
Deltav is the peculiar line-of-sight velocity: vgal-vc
vc is cluster velocity 
"""

#Others
G = 6.6738480E-11
Mpc2m = 3.0857E+22
M_sun = 1.98847E30 #kg

#Define parameters from 11_HydraRSP
ind = 0
table = Table()
table['sigma'][ind] = 620    #km/s
table['r200'][ind] = 1.4   #Mpc
table['d_HC'][ind] = 58.6   #Mpc
table['v_c'][ind] = 3686   #km/s

table['M200'][ind] = 3.02E14*M_sun

table['ta'][ind] = 0.5          #kpc
table['Rc'][ind] = 65             #kpc
table['Rc_Hydra'][ind] = 15.94    #kpc
table['rho0'][ind]=14.7E-3        #cm^-3
table['rho0_Hydra'][ind]=1.49E-25 #g/cm^3
table['Md_s'][ind]=M_sun*4.6E10   #kg
table['Md_g'][ind]=0.06*Md_s      #kg
table['Rd_s'][ind]=2.15           #kpc
table['Rd_g'][ind]=1.7*Rd_s       #kpc
table['kpc2m'][ind] = 3.0857E+19  #m          #There are 2 magn order missing if I do not use E+20. Maybe is Rd.
table['sigmas'][ind] = 620000      #m/s
table['edge'][ind] = 2.7

"""
Function that estimates the stripping velocities from Wang21

Libraries:  import astropy.units as u
            import numpy as np
            from astropy.constants import c
            from astropy.table import Table

Inputs: table -> astropy table with the parameters of the cluster(s)
        ind -> row of the cluster inside table
        
Outputs: vstripx -> absolute value of stripping velocity

To improve -> Add units to tables to get rid of conversion factors
"""

def f_vstrip(table,ind):
    #-------------------------------------------------------
    Sigma_0s = table['Md_s'][ind]/(2*np.pi*table['Rd_s'][ind]*table['Rd_s'][ind]*kpc2m*kpc2m)
    Sigma_s = Sigma_0s*np.exp(-6.1/table['Rd_s'][ind])
    Sigma_0g = table['Md_g'][ind]/(2*np.pi*table['Rd_g'][ind]*table['Rd_g'][ind]*kpc2m*kpc2m)
    Sigma_g = Sigma_0g*np.exp(-6.1/table['Rd_g'][ind])
    Pi = 2*np.pi*G*table['Sigma_s'][ind]*table['Sigma_g'][ind]
    #------------------------------------------------------
    xs = np.linspace(0.002,table['edge'][ind])
    rho_icm = table['rho0_Hydra'][ind]*1.0E3*(1+(0.5*np.pi*xs*1000*table['r200'][ind]*1000/table['Rc'][ind])**2)**(-3*table['beta'][ind]/2)
    ys = np.sqrt(Pi/(3*rho_icm))
    #-------------------------------------------------------
    vstrip = ys/table['sigmas']

    return(vstrip)

"""
Function that estimates the escape velocities from Wang21

Libraries:  import astropy.units as u
            import numpy as np
            from astropy.constants import c

Inputs: params -> astropy table with the parameters of the cluster(s)
        
Outputs: vesc -> absolute value of scape velocity
To improve -> Add units to tables to get rid of conversion factors
"""
def f_vesc(table,ind):
    #Define vesc
    x1 = np.linspace(0.002,table['edge'][ind]/2)
    x2 = np.linspace(table['edge'][ind]/2,table['edge'][ind])
    #---------------------------------------
    c = 6
    s1 = 0.5*np.pi*x1/table['r200'][ind]
    s2 = 0.5*np.pi*x2/table['r200'][ind]
    gc = 1/(np.log(1+c)-c/(1+c))
    K1 = gc*(np.log(1+c*s1)/s1-np.log(1+c))+1
    K2 = gc*(np.log(1+c*s2)/s2-np.log(1+c))+1
    y1 = np.sqrt(2*G*table['M200'][ind]*K1/(3*table['r200'][ind]*Mpc2m))
    y2 = np.sqrt(2*G*table['M200'][ind]/(3*table['r200'][ind]*s2*Mpc2m))
    #----------------------------------------
    dif=y2[0]-y1[49]
    y2 = y2-dif
    x = np.concatenate((x1,x2))
    y = np.concatenate((y1,y2))
    #-----------------------------------------
    vesc = y/(table['sigma'][ind]*1E3)

    return(vesc)

