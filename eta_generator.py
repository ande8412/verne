import sys
sys.path.append('./src/')
from src.CalcEtaDist_light import *
import numpy as np
from tqdm.autonotebook import tqdm
# heavy_sigma_list = np.logspace(-42,-32,num=20)
# flipped_heavy_sigma_list = np.flip(heavy_sigma_list)
# heavy_dm_masses = np.concatenate((np.array([1,5,10,20]),np.arange(20,50,5),np.arange(50,110,10),np.arange(200,1000,20),np.array([1000,1500,2000])))


# sigma_list =np.logspace(-40,-27,25)
# flipped_sigma_list = np.flip(sigma_list)
# dm_masses = np.concatenate((np.arange(0.1,2,0.1),np.arange(2,5,0.5),np.arange(5,11,0.5),np.array([20,50,100,1000])))

# fine_sigma_list =np.logspace(-35,-29,50)
# flipped_fine_sigma_list =np.flip(fine_sigma_list)
# fine_dm_masses = np.concatenate((np.arange(0.6,2,0.05),np.arange(2,5,0.1)))#,np.arange(5,11,0.5)))#,np.array([20,50,100,1000])))


# heavy_sigma_list2 = np.logspace(-39,-32,num=20)
# heavy_dm_masses2 = np.concatenate((np.arange(20,50,5),np.arange(50,110,5),np.arange(200,1000,10),np.array([1000,1500,2000])))
# flipped_heavy_sigma_list2 = np.flip(heavy_sigma_list2)
# gen_eta(flipped_heavy_sigma_list2,heavy_dm_masses2,overwrite=False)
sigmas = np.logspace(-43,-29,num=50)
light_masses = np.geomspace(0.1,1,50)
medium_masses = np.geomspace(1,100,50)
flipped_sigmas = np.flip(sigmas)
heavy_masses = np.geomspace(100,1000,50)


mX_grid = np.geomspace(0.5,1000,50)
mX_grid = np.concatenate((mX_grid,np.array([0.6,1,10,100,1000])))
mX_grid = np.unique(np.sort(mX_grid))

sigma_grid = np.array([1e-40,1e-39,1e-38,1e-37,1e-36,1e-35])
flipped_sigma_grid = np.flip(sigma_grid)


heavy_mediator_heavymX_sigmaE_grid = [1e-42,1e-41,1e-40,1e-39,1e-38,1e-37,1e-36]
heavy_mediator_lightmX_sigmaE_grid =  [1e-39,3.16e-39,1e-38,3.16e-38,1e-37,3.16e-37,1e-36,3.16e-36,1e-35,3.16e-35,1e-34,3.16e-34,1e-33,3.16e-33,1e-32,3.16e-32,1e-31,3.16e-31,1e-30,3.16e-30,1e-29]

light_mediator_heavymX_sigmaE_grid = [1e-37,1e-36,1e-35,1e-34,1e-33]
light_mediator_lightmX_sigmaE_grid =  [1e-37,3.16e-37,1e-36,3.16e-36,1e-35,3.16e-35,1e-34,3.16e-34,1e-33,3.16e-33,1e-32,3.16e-32,1e-31,3.16e-31,1e-30,3.16e-30,1e-29]


heavymX_grid = np.geomspace(10,100,50)
heavymX_grid = np.concatenate((heavymX_grid,np.array([10,100,1000])))
heavymX_grid = np.unique(np.sort(heavymX_grid))

lightmX_grid = np.geomspace(0.5,10,50)
lightmX_grid = np.concatenate((lightmX_grid,np.array([0.6,1,10])))
lightmX_grid = np.unique(np.sort(lightmX_grid))

lightmX_damascus_matching = np.geomspace(0.5,10,25)
lightmX_damascus_matching = np.concatenate((lightmX_damascus_matching,np.array([0.6,1,10])))
lightmX_damascus_matching = np.unique(np.sort(lightmX_damascus_matching))


heavier_mX_grid = np.geomspace(100,1000,50)
# heavier_mX_grid = np.unique(np.sort(heavier_mX_grid))


full_mX_grid = np.geomspace(0.5,1000,50)
full_mX_grid = np.concatenate((full_mX_grid,np.array([0.6,1,10,20,30,40,50,100,1000])))
full_mX_grid = np.unique(np.sort(full_mX_grid))

full_sE_grid = [1e-44,3.16e-44,1e-43,3.16e-43,1e-42,3.16e-42,1e-41,3.16e-41,1e-40,3.16e-40,1e-39,3.16e-39,1e-38,3.16e-38,1e-37,3.16e-37,1e-36,3.16e-36,1e-35,3.16e-35,1e-34,3.16e-34,1e-33,3.16e-33,1e-32,3.16e-32,1e-31,3.16e-31,1e-30,3.16e-30,1e-29]

def gen_eta(sigma_list,dm_masses,overwrite=False):
    qedir = '/Users/ansh/Local/SENSEI/sensei_toy_limit/python/theory_tools/QEDark/halo_data/modulated/'
    for mediator in ['hm','ulm']:
        for s in tqdm(range(len(sigma_list))):
            for m in tqdm(range(len(dm_masses))):
                sigmaE = sigma_list[s]
                sigmaE = float(format(sigmaE, '.3g'))
                mX = dm_masses[m]
                mX = np.round(mX,2)
                calcEtaDist_light(mX, sigmaE, loc='full', interaction=mediator, depth = 104,num_angles=36,outdir=qedir,overwrite=overwrite,write=True)


#  gen_eta(sigmas,light_masses,overwrite=False)
#  gen_eta(flipped_sigmas,light_masses,overwrite=False)
#  gen_eta(sigmas,heavy_masses,overwrite=False)
#  gen_eta(flipped_sigmas,heavy_masses,overwrite=False)
#  gen_eta(sigmas,medium_masses,overwrite=False)
#  gen_eta(flipped_sigmas,medium_masses,overwrite=False)
# gen_eta(heavy_mediator_heavymX_sigmaE_grid,heavymX_grid,overwrite=False)
# gen_eta(light_mediator_heavymX_sigmaE_grid,heavymX_grid,overwrite=False)


# gen_eta(heavy_mediator_lightmX_sigmaE_grid,lightmX_grid,overwrite=False)
# gen_eta(light_mediator_lightmX_sigmaE_grid,lightmX_grid,overwrite=False)

# gen_eta(heavy_mediator_lightmX_sigmaE_grid,lightmX_damascus_matching,overwrite=False)
# gen_eta(light_mediator_lightmX_sigmaE_grid,lightmX_damascus_matching,overwrite=False)

# gen_eta(heavy_mediator_heavymX_sigmaE_grid,heavier_mX_grid,overwrite=False)

# gen_eta(full_sE_grid,full_mX_grid,overwrite=False)

# gen_eta(light_mediator_heavymX_sigmaE_grid,heavier_mX_grid,overwrite=False)

#sigmalist = np.geomspace(1e-29,1e-35,20)