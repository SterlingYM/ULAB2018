import numpy as np
import pylab as pl
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
fname = 'LFI_SkyMap_030_1024_R2.01_full.fits'
tmap= hp.read_map(fname)
rawdata= fits.open(fname)
data=hp.visufunc.mollview(tmap, return_projected_map= True)
print(rawdata)
