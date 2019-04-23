import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy_healpix as ap_h

cmap = 'rainbow'
NSIDE = 1024
fname = "LFI_SkyMap_030_1024_R2.01_full.fits"

# index to coord
hp_index = np.arange(hp.nside2npix(NSIDE))
lon,lat = ap_h.healpix_to_lonlat(hp_index,nside=NSIDE,order='ring')

# get data
healpix_data = hp.read_map(fname)

# Galactic coord data:
x = lon.degree
for i in range(len(x)):
    if x[i]>180:
        x[i] = x[i]-360
x = np.radians(x)
y = np.radians(lat.degree)
z = healpix_data

galactic_data = [x,y,z] # <--------data[lon(radian),lat(radian),value)]


# Ecliptic coord data:
#for some reason this rotator is not working properly!
x_e,y_e = hp.Rotator(coord='ge',deg=True)(np.degrees(x),np.degrees(y)) # rotation
print(x_e.min())
print(x_e.max())
print(y_e.min())
print(y_e.max())
for i in range(len(y_e)):
    y_e[i] += np.pi/2
    if y_e[i]>np.pi/2:
        y_e[i]=  y_e[i] - np.pi
z_e = healpix_data
ecliptic_data = [x_e,y_e,z_e] # <--------data


# plot test
plt.figure(figsize=(15,10))
ax = plt.subplot(111,projection='mollweide')
ax.scatter(x[0::30],y[0::30],c=np.log(z[0::30])*100,s=10)
plt.draw()

plt.figure(figsize=(15,10))
ax = plt.subplot(111,projection='mollweide')
ax.scatter(x_e[0::30],y_e[0::30],c=np.log(z_e[0::30])*100,s=10)
plt.draw()



#
#ecliptic_healpix_data = hp.Rotator(coord=['G','E'])
#
#
## in original frame
#hp.visufunc.mollview(
#    data,
#    title="Histogram equalized Galactic",
#    norm="hist",
#    cmap=cmap,
#    return_projected_map=True
#)
#plt.draw()
#
## in ecliptic coord.
hp.visufunc.mollview(
    healpix_data,
    coord=["G", "E"],
    title="Histogram equalized Ecliptic",
    norm="hist",
    cmap=cmap,
    return_projected_map=True
)
plt.draw()


plt.show()
