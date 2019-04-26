# CMB_sampling.py
#
# Takes Planck CMB data fits file and converts into array in various ways:
#   (1) galactic coord with max sample (get_galactic())
#   (2) ecliptic coord with max sample (get_ecliptic())
#   (3) at specified resolution (sampling())
#   (4) scatter test plot (plot_test())
#   (5) rectangular plot(plot_rect()) # complete version of this function is in CMB_vs_H0.py
#
# author: Yukei Murakami (sterling.astro@berkeley.edu)
# last edit: 4/26/2019

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy_healpix as ap_h
import pickle
from multiprocessing import Pool
import multiprocessing
from functools import partial

multiprocess = True # change to turn on/off multicore processing
N_SAMPLES = 60
NSIDE = 1024
fname = "LFI_SkyMap_030_1024_R2.01_full.fits"
save_filename = "CMB_sampled.dat"

def get_lonlat(NSIDE):
    # index to coord
    hp_index = np.arange(hp.nside2npix(NSIDE))
    lon,lat = ap_h.healpix_to_lonlat(hp_index,nside=NSIDE,order='ring')
    x = lon.degree
    for i in range(len(x)):
        if x[i]>180:
            x[i] = x[i]-360
    x = np.radians(x)
    y = np.radians(lat.degree)
    return x,y

def get_galactic(fname,NSIDE):
    # Galactic coord data:
    # returns array of [[radians,radians,value], . . . ]
    print('Generating data (galactic): it may take a while.')
    x,y = get_lonlat(NSIDE)
    z = hp.read_map(fname)
    return [x,y,z]

def get_ecliptic(fname,NSIDE):
    # Ecliptic coord data:
    # returns array of [[radians,radians,value], . . . ]
    print('Generating data (ecliptic): it may take a while.')
    x,y = get_lonlat(NSIDE)
    z_e = hp.Rotator(coord='ge',deg=False).rotate_map_alms(hp.read_map(fname))
    return [x,y,z_e]

def x_search(x,y,z,x_low,x_i,y_low,y_i):
    local_total = 0
    local_count = 0
    for k in range(len(z)):
        if (x_low <= x[k] and x[k] < (x_low + x_i))\
                and (y_low <= y[k] and y[k] < (y_low + y_i)):
            local_total += z[k]
            local_count += 1
    if local_count != 0:
        local_avg = local_total / local_count
    else:
        local_avg = 0
    local_data = [x_low,x_low+x_i,y_low,y_low+y_i,local_avg]
    return local_data

def sampling(data,n_rows,n_cols):
    # data
    x = data[0]
    y = data[1]
    z = data[2]
    
    # prep
    x_i = 2*np.pi/n_cols # x increment
    y_i = np.pi/n_rows # y increment
    x_low = -np.pi # initial lower bound
    y_low = -np.pi/2 # initial lower bound
     
    # loop
    sampled_data = []
    for i in range(n_rows):
        x_low = -np.pi # initial lower bound
        x_low_list = np.arange(x_low,x_i*n_cols,x_i)   
        
        if multiprocess == True:
            with Pool(7) as p:
                partial_map = partial(x_search, x, y, z, x_i=x_i, y_low=y_low, y_i=y_i)
                x_search_data = p.map(partial_map, x_low_list)
                sampled_data.extend(x_search_data)
        else:
            for x_low in x_low_list:
                local_total = 0
                local_count = 0
                for k in range(len(z)):
                    if (x_low <= x[k] and x[k] < (x_low + x_i))\
                            and (y_low <= y[k] and y[k] < (y_low + y_i)):
                        local_total += z[k]
                        local_count += 1
                if local_count != 0:
                    local_avg = local_total / local_count
                else:
                    local_avg = 0
                local_data = [x_low,x_low+x_i,y_low,y_low+y_i,local_avg]
                sampled_data.append(local_data)
 
        y_low += y_i
        print('\rSampling Data: {:.1f}% done '.format(i*100/n_rows),end='')
    print('')
	
    # return data: [[rad,rad,rad,rad,val], . . . ]
    return sampled_data

def plot_test(data,niter=30):
    # plot test: plots in mollweide projection with log scaling
    # increase niter for faster plot (plots every 'niter'th point
    # (often large data slows down plotting)
    plt.figure(figsize=(10,6))
    ax = plt.subplot(111,projection='mollweide')
    ax.scatter(data[0][0::niter],data[1][0::niter],c=np.log10(data[2][0::niter]),s=1)
    plt.draw()

def plot_rect(data):
    data = np.array(data)
    data = data.transpose()

    #print(data)
    plt.figure(figsize=(10,6))
    ax = plt.subplot(111,projection='mollweide')
    ax.scatter((data[0]+data[1])/2,(data[2]+data[3])/2,c=np.log10(data[4]),s=100)
    plt.draw()

def savedata(save_filename,data):
    output=open(save_filename,'wb')
    pickle.dump(data,output)
    output.close()

# main
#plot_test(get_galactic(fname,NSIDE),niter=30)
#plot_test(get_ecliptic(fname,NSIDE),niter=30)

sampled_data = sampling(get_ecliptic(fname,NSIDE),N_SAMPLES,N_SAMPLES)
savedata(save_filename,sampled_data)
plot_rect(sampled_data)

plt.show()
