import numpy as np
import healpy as hp
import math
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

alphaG= math.radians(192.86)
deltaG=math.radians(27.13)
LNCP=math.radians(122.93)

def gal_to_ecl(lon,lat):
    declination=np.arcsin(math.sin(deltaG)*math.sin(lat)+math.cos(deltaG)*math.cos(lat)*math.cos(LNCP-lon))
    acension=     np.arcsin( math.sin(LNCP-lon)*math.cos(lat)/math.cos(declination))+alphaG
    return [declination,acension]






temp=healpix_data

declination=[]
acension=[]

for i in range(len(x)):
    declination.append(gal_to_ecl(x[i],y[i])[0])
    acension.append(gal_to_ecl(x[i],y[i])[1])


# Ecliptic coord data:
#for some reason this rotator is not working properly!
print(min(x_e))
print(max(x_e))
print(min(y_e))
print(max(y_e))
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
ax.scatter(declination,acension,c=np.log(temp[0::30])*100,s=10)
plt.draw()



#
ecliptic_healpix_data = hp.Rotator(coord=['G','E'])
#
#
## in original frame
hp.visufunc.mollview(
    data,
    title="Histogram equalized Galactic",
    norm="hist",
    cmap=cmap,
    return_projected_map=True
)
plt.draw()
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

plt.show()
