import numpy as np
import pylab as pl
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
fname = 'LFI_SkyMap_030_1024_R2.01_full.fits'
tmap= hp.read_map(fname)
hp.visufunc.mollview(tmap)
plt.show()
hdulist= fits.open(fname)



rawdata=np.array(hdulist[1].data)
print(rawdata)
