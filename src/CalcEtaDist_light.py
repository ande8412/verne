import numpy as np
from scipy.interpolate import interp1d

import MaxwellBoltzmann as MB
import argparse
import os
import verne
import verne_light

from matplotlib import pyplot as plt

try:
    from tqdm import tqdm
except:
    def tqdm(x):
        return x
    
me_eV = 5.1099894e5
mP_eV = 938.27208816 *1e6
def mu_Xe(mX):
    """
    DM-electron reduced mass
    """
    return mX*me_eV/(mX+me_eV)


def mu_XP(mX):
    """
    DM-proton reduced mass
    """
    return mX*mP_eV/(mX+mP_eV)

def sigmaE_to_sigmaP(sigmaE,mX):
    import numpy as np
    mX*=1e6 #eV
    sigmaP = sigmaE*(mu_XP(mX)/mu_Xe(mX))**2
    # sigmaP = np.rou
    return sigmaP


def etaIntegral(v, interpfun):
        from scipy.integrate import quad
        return quad(lambda x: interpfun(x)/x, v, 800.0, epsrel=1e-3)[0]


def calcEtaDist_light(mX, sigmaE, loc, interaction, depth = 0,num_angles=36,outdir='./',overwrite=False,write=True):
    
    
    import csv
    
    if (interaction not in  ["SI", "hm", "ulm"]):
        print("> Unknown interaction type <", interaction, ">...")
        exit()

    if interaction == 'hm':
        Mediator = 'Scr'
    elif interaction == 'ulm':
        Mediator = 'LM'
    else:
        print('Not applicable for DM Electron Scattering')
         
    target = loc
    if (target.lower() not in ["atmos", "full", "earth"]):
        print("> Only allowed targets for light DM are 'atmos', 'full', 'earth'...")
        exit()


    


    mX_str = '{0:.4f}MeV'.format(mX) 
    sigmaE = float(format(sigmaE, '.3g'))
    sigmaP = sigmaE_to_sigmaP(sigmaE,mX)
    sigmaP = float(format(sigmaP, '.3g'))

    m_x = mX *1e-3
    sigma_p = sigmaP
    mX = float(mX)
    mass_str = str(np.round(mX,3)).replace('.','_')
    print(" ")
    print("> Calculating for...")
    print(">    m_x/GeV:", m_x)
    print(">    sigma_p/cm^2:", sigma_p)
    #print "        gamma/pi:", gamma_by_pi
    print(">    detector at :", loc)
    print(" ")


    if write:
        write_dir  = outdir + f'Verne_{Mediator}/mDM_{mass_str}_MeV_sigmaE_{sigmaE}_cm2/'
        if not os.path.isdir(write_dir):
            os.mkdir(write_dir)

        num_angles_present = len(os.listdir(write_dir))
        if num_angles_present != num_angles and num_angles_present > 0:
            print(f'num isoangles in {write_dir} dont match the one set by this code! Overwriting')
            overwrite = True
        if overwrite:
            print(f'Overwrite Flag is on, deleting all angles in: {write_dir}')
            files = os.listdir(write_dir)
            for file in files:
                print(file)
                os.remove(write_dir + file)
        if not overwrite and num_angles_present == num_angles:
            print("point generated, continuing")
            return 0
    
    
    #Initialise verne
    verne_light.Tabulate_Column_Density(depth, target)
            
    #Loop over gamma values
    N_gamma = num_angles
    Nv1 = 49 #Near the velocity threshold
    Nv2 = 50 #Everywhere else
    Nv = Nv1 + Nv2  + 1 + 1 #Add an extra one for 20 km/s
            
    v_th = 20e0 #Lowest speed to consider
    vmax = MB.ve + MB.vesc
    
    def getVelDist(gamma):
        
        #Generate a list of sampling values for v (with some very close to v_th)
        #vlist = np.geomspace(v_th, 0.25*vmax, Nv1)    
        #vlist = np.append(vlist, np.linspace(0.15*vmax, 0.89*vmax, Nv2)) 
        #vlist = np.append(vlist, np.linspace(0.9*vmax, 0.9999*vmax, Nv3)) 
        
        vlist = np.linspace(v_th, 0.84*vmax, Nv1)
        vlist = np.append(vlist, np.linspace(0.85, 1.0, Nv2)*vmax)
        vlist = np.append(vlist, 0.99*v_th)
        vlist = np.sort(vlist)
        
        f_final = 0.0*vlist
        for i in range(len(vlist)):
            f_trans = verne_light.CalcF_transmitted(vlist[i], gamma, sigma_p, m_x, target, interaction)
            f_refl = verne_light.CalcF_reflected(vlist[i], gamma, sigma_p, m_x, target, interaction)
            f_final[i] = f_trans + f_refl
    
        #Add on the final point
        vlist = np.append(vlist, vmax)
        f_final = np.append(f_final, 0.0)
    
        return vlist, f_final
    
    gamma_list = np.linspace(0, 1, N_gamma)
    gamma_list[0] = 1e-3
    gamma_list[-1] = 1 - 1e-3

    gamma_rep = np.repeat(gamma_list, Nv)

    vgrid = np.zeros((N_gamma, Nv))
    fgrid = np.zeros((N_gamma, Nv))
    fgrid_withrefl = np.zeros((N_gamma, Nv))
    etas_by_angle = []
    vs_by_angle = []
    fvals_by_angle = []

    for j in tqdm(range(N_gamma), desc='Calculating velocity distribution'):
        etas = []
        #print(">Calculating for gamma/pi = ", gamma_list[j],"...")
        vvals,fvals = getVelDist(gamma_list[j]*np.pi)
        # vgrid[j,:], fgrid[j,:] =  vvals,fvals
        f_interp = interp1d(vvals, fvals, kind='linear',bounds_error=False, fill_value=0.0)
        vs_by_angle.append(vvals)
        for v in vvals:
            eta = etaIntegral(v, f_interp)
            etas.append(eta)
        etas_by_angle.append(np.array(etas))
        fvals_by_angle.append(np.array(fvals))
        
    etagrid = np.array(etas_by_angle)
    vgrid = np.array(vs_by_angle)
    fvgrid = np.array(fvals_by_angle)
    


    
    
    
    if write:


        isoangles = np.arange(num_angles-1,-1,-1)
        # print('writing')
        for i in range(num_angles):
            isoangle = isoangles[i]
            # print('isoangle:',isoangle)
            with open(write_dir+f'DM_Eta_theta_{i}.txt','w') as f:
                writer = csv.writer(f,delimiter='\t')
                v = vgrid[isoangle,:]
                eta = etagrid[isoangle,:]
                writer.writerows(zip(v,eta))

        
    else:
        return vgrid,etagrid,fvgrid