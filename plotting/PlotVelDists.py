import sys
sys.path.append("../src")

import numpy as np
import MaxwellBoltzmann as MB

#Matplotlib ------------
import matplotlib
font = { 'size'   : 16, 'family':'serif'}
matplotlib.rc('font', **font)

import matplotlib.pyplot as plt
#------------------------

def MakePlot(filename, N_gamma, N_v, logy = False):
    
    #Load velocity distributions and reshape
    _gammas, _vs, _fs = np.loadtxt(filename, unpack=True)
    gamma_grid = np.reshape(_gammas, (N_gamma, N_v))
    v_grid = np.reshape(_vs, (N_gamma, N_v))
    f_grid = np.reshape(_fs, (N_gamma, N_v))
    

    
    colmap = matplotlib.cm.get_cmap('turbo_r')
    
    plt.figure()
    
    #Plot the 'unperturbed' Maxwell-Boltzmann Distribution
    v_list = np.linspace(0.0, 800.0, 1000)
    plt.plot(v_list, MB.calcf_SHM(v_list), linestyle='--', color='k')

    #Plot velocity distributions for different values of gamma = pi - Theta
    for i in range(N_gamma):
        plt.plot(v_grid[i,:], f_grid[i,:], color=colmap(i/(N_gamma-1)), alpha = 0.5, linestyle='-')


    if (logy):
        plt.ylim(1e-5, 1e-2)
        plt.yscale('log')
    else:
        plt.ylim(0, 0.005)

    plt.xlim(0, 800) 

    plt.xlabel(r'$v$ [km/s]')
    plt.ylabel(r'$f_\chi(v)$ [(km/s)$^{-1}$]')

    plt.plot(-100, 0, color=colmap(0.99), label=r'$\Theta = 0^\circ$')
    plt.plot(-100, 0, color=colmap(0.75), label=r'$\Theta = 45^\circ$')
    plt.plot(-100, 0, color=colmap(0.5), label=r'$\Theta = 90^\circ$')
    plt.plot(-100, 0, color=colmap(0.25), label=r'$\Theta = 135^\circ$')
    plt.plot(-100, 0, color=colmap(0.01), label=r'$\Theta = 180^\circ$')

    plt.legend(loc='best')

    plt.savefig("../plots/VelDist.pdf")
    plt.show()
    
filename = "f_light_hm_full_mx10.0000MeV_lsig-34.38.txt"
MakePlot("../results/veldists/" + filename, N_gamma=41, N_v=100, logy = False)