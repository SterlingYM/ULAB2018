import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.optimize as fitter
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap
import pandas as pd
#import healpy as hp

##
#@source: 1.Open Supernova Catalog
#   2.https://www.cfa.harvard.edu/supernova/README_archive
#   3.http://balbuceosastropy.blogspot.com/2013/09/the-mollweide-projection.html
#   4.https://zonca.github.io/2013/09/Planck-CMB-map-at-high-resolution.html
#@author: Xuyang Yu
#plot the distribution of Type Ia supernova on a XY-grid,
#where x-axis is the right ascension and y-axis is the declination
#domain: x is [0, 360], y is [-90, 90]

#Ia-norm    : spectroscopically/photometrically "normal" type-Ia supernovae
# Ia-91t     : type-Ia supernovae whose lightcurves/spectra resemble those of SN 1991T
# Ia-91bg    : type-Ia supernovae whose lightcurves/spectra resemble those of SN 1991bg
# Ia-pec     : all other type-Ia supernovae (e.g., 2000cx)
#
#did not add in the correction from 3.2E(B-V) but not sure if it is right
#add in a filter of years >= 2005
#disabled the repeating filter (True or...)
#disabled the name Ia filter
#Linear fit: Hubble constants: 6.329e-05 data points: 8651
#contrast enhancement centered around 0.036 +- 0.004



# def hasColumn(row, index):
#     if row[index] == '':
#         return False
#     return True

def hasAllColumns(row, totColNum):
    for i in range(totColNum):
        if row[i] == '':
            return False
    #if row[SNeTypeIND] != 'Ia':
        #return False
    year = (int) (row[dateIND][0:4])
    if year < 2005:
        return False
    return True

def readIn(fileName):
    name = []
    date = []
    app_mag = []
    abs_mag = []
    ra = []
    dec = []
    z = []
    SNeType = []
    ebv = []
    csv_file = open(fileName, 'r')
    data_reader = csv.reader(csv_file, delimiter = ',')
    next(data_reader)
    for row in data_reader:
        if hasAllColumns(row, totColNum):
            name.append(row[nameIND])
            date.append(row[dateIND])
            app_mag.append(row[app_magIND])
            abs_mag.append(row[abs_magIND])
            ra.append(row[raIND])
            dec.append(row[decIND])
            z.append(row[zIND])
            SNeType.append(row[SNeTypeIND])
            ebv.append(row[ebvIND])
    return name, date, app_mag, abs_mag, ra, dec, z, SNeType, ebv

def str2float(data_list):
    for i in range(len(data_list)):
        data_list[i] = float(data_list[i])
    return data_list

def dist_velocity(name, app_mag, abs_mag, z, ebv):
    dist = []
    dist_rej = []
    v = []
    v_rej = []
    approved_SNe_data = []
    rejected_SNe_data = []

    for i in range(len(z)):
        x = 10 ** ((app_mag[i] - abs_mag[i] - 0 * ebv[i] + 5) / 5)
        y = z[i] * c
        if True or y not in v:
            dist.append(x)
            v.append(y)
            approved_SNe_data.append([name[i], app_mag[i], abs_mag[i], z[i]])
        else:
            rejected_SNe_data.append([name[i],app_mag[i],abs_mag[i],z[i]])
            dist_rej.append(x)
            v_rej.append(y)
    return approved_SNe_data, dist, v, rejected_SNe_data, dist_rej, v_rej



c = 3.0 * 10 ** 5 #unit [km/s]

#if it does not exist, put -1; if multiple value, take arithmetic mean
def purifyValues(v, isDeg, hasSign):
    for i in range(len(v)):
        if v[i] == '':
           v[i] = -1
        else:
            v_a = 0.0
            strList = v[i].split(',')
            for j in range(len(strList)):
                if isDeg:
                    v_a = v_a + degMinSecToDeg(strList[j], hasSign)
                else:
                    v_a = v_a + (float)(strList[j])
            v[i] = v_a / len(strList)
    return v

def degMinSecToDeg(str, hasSign):
    strList = str.split(':')
    multiplier = 1
    if hasSign:
        sign = (strList[0])[:1]
        strList[0] = (strList[0])[1:]
        if sign == '-':
            multiplier = -1
    deg = 0.0
    for i in range(len(strList)):
        deg += (float)(strList[i]) / (60 ** i)
    return deg * multiplier

def raToDegree(ra):
    for i in range(len(ra)):
        ra[i] = ra[i] * 15
    return ra

#pls make sure 1st row of data is valid!
#return max and min of non-negatives
def findMaxMin(v):
    min = v[0]
    max = v[0]
    for i in range(len(v)):
        if v[i] < 0:
            continue
        if min > v[i]:
            min = v[i]
        if max < v[i]:
            max = v[i]
    return max, min

def findZeroDec(dec):
    count = 0
    for i in range(len(dec)):
        if dec[i] < zeroDecUp and dec[i] > zeroDecDown:
            count = count + 1
    return count

fileName = 'Name2E(B-V).csv'

# def calculateXY(ra, dec):
#     x = []
#     y = []
#     r = 100
#     for i in range(len(ra)):
#         phi = dec[i]*np.pi/180
#         theta = NewtonRaphsonIteration(phi, phi, 30)
#         xVal = r*2*np.sqrt(2)*(ra[i]-180)/180*np.cos(theta)
#         yVal = r*np.sqrt(2)*np.sin(theta)
#         x.append(xVal)
#         y.append(yVal)
#     return x, y

def degreeToRadian(ra, dec):
    raRad = []
    decRad = []
    for i in range(len(ra)):
        raRad.append(ra[i]*np.pi/180)
    for i in range(len(dec)):
        decRad.append((dec[i] - 180) * np.pi/180)
    return raRad, decRad


def plotLinearH():
    approved_data, dist, v, rejected_data, dist_rej, v_rej = dist_velocity(name, app_mag, abs_mag, z, ebv)
    print("Total number of data: " + str(len(v)))
    print("Total number of rejected data: " + str(len(v_rej)))
    #plt.figure(figsize = (15,5))
    #plt.scatter(dist, v, s = 5, c = 'blue')
    #plt.show()

    print("\n\nPart 2: fitting")
    ave_z = fitting(my_model0, dist, v)
    return ave_z, dist, v

def createPoints(ra, dec):
    points = []
    for i in range(len(ra)):
        points.append([ra[i], dec[i]])
    return points

# def interpolateH(ra, dec, Color): #smooth interpolation
#     points = createPoints(ra, dec)
#     grid_x, grid_y = np.mgrid[0:1:360j, 0:1:180j]
#     grid = griddata(points, Color, (grid_x, grid_y), method='cubic')
#     plt.figure()
#     plt.imshow(grid, extent=(0,1,0,1), origin='lower', interpolation= 'bicubic')
#     np.savetxt("interpolatedGrid.cvs", grid, delimiter=",")
#     totIntDataNum = 360 * 180
#     newRA = []
#     newDEC = []
#     for i in range(totIntDataNum):
#         newRA.append(i%360)
#         newDEC.append(-(i//360) + 90)
#     newColor = []
#     for row in grid:
#         for val in row:
#             if val == 'nan':
#                 newColor.append(0.04)
#             else:
#                 newColor.append(val)
#     plot_mwd(newRA, newDEC, newColor)
#     plt.show()

#upFilter = 0.043
#downFilter = 0.035
upFilter = 0.1
downFilter = -0.1

def calculateColor(ave_z, dist, v):
    Color = []
    for i in range((len(dist))):
        Color.append((v[i]/ dist[i] - ave_z) / ave_z)

    max = 0
    min = 1000000000
    for i in range(len(Color)):
        if Color[i] > max:
            max = Color[i]
        if Color[i] < min:
            min = Color[i]
        #(after normalization)
    if False:
        for i in range(len(Color)):
            Color[i] = (Color[i] - min) / (max - min)
    over = 0
    under = 0
    if True:
        for i in range(len(Color)):
            if Color[i] > upFilter:
                Color[i] = upFilter
                over = over + 1
            if Color[i] < downFilter:
                Color[i] = downFilter
                under = under + 1
    print(max)
    print(len(Color))
    print(over)
    print(under)
    for i in range(4):
        print(Color[i])

    np.savetxt("calculatedColor.cvs", Color, delimiter=",")

    return Color

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
       return v
    return v / norm

def my_model0(x, H):
    return H * x

# def my_model1(x, H, k):
#     return H * x + k

# def my_model2(x, H, a):
#     return H * x + a * x ** 2

def fitting(model, xdata, ydata):
    # data fitting
    par0 = np.array([0.01]) # initial guess
    par, cov = fitter.curve_fit(model, xdata, ydata, par0, absolute_sigma=True)

    # plot
    x = np.linspace(0, max(xdata), 10**4)
    y_fitted = my_model0(x, par[0])

    if False:
        plt.figure(figsize = (15,5))
        plt.plot(x,y_fitted, c='orange')
        plt.scatter(xdata, ydata, s=5)
        plt.title("Fitted function")
        plt.legend(['da/dt = {:2f}a'.format(par[0])])
        plt.xlabel('distance [pc]')
        plt.ylabel('velocity [km/s]')
    #plt.show()

    # result
    print("Hubble Constant H0 = {:.4} [km/(s*pc)]".format(par[0]))
    return par[0]
    #print("y-intercept: " + str(par[1]))
    #print("accel: " + str(par[1]))

def plot_mwd(RA, Dec, Color, ifFillRect, org=0, title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = [np.remainder(n+360-org,360) for n in RA] # shift RA values
    for i in range(len(x)):
        if x[i] > 180:
            x[i] -=360    # scale conversion to [-180, 180]
        x[i] =-x[i]    # reverse the scale: East to the left
    colombi1_cmap = ListedColormap(np.loadtxt("CMBColorMap.txt")/255.)
    colombi1_cmap.set_bad("gray") # color of missing pixels
    colombi1_cmap.set_under("white") # color of background, necessary if you want to use
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection, facecolor ='darkgray')
    rand = np.random.random_sample((8651,))
    ax.scatter(np.radians(x), np.radians(Dec), c = Color, s = 1, alpha=1, cmap= colombi1_cmap)  # convert degrees to radians
    #fig.colorbar(ax, orientation='horizontal', fraction=.1)
    #ax.scatter(0, 0, c='red', s = 10)
    if ifFillRect:
        fillRect(x, Dec, Color, 6, 3, ax)

    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
    plt.show()



def fillRect(ra, dec, Color, rectW, rectH, ax): #sliced average
    whiteOutLimit = 1 #if num < whiteOutLimit, the region becomes white
    fillGapW = rectW/2000
    fillGapH = rectH/100
    totCol = (int)(360/rectW)
    totRow = (int)(180/rectH)
    colorArray = np.loadtxt("CMBColorMap.txt")/255
    print("rectangleNum: " + (str)(totCol*totRow))
    rectColor = [x[:] for x in [[0] * totCol] * totRow]
    rectColorPointCount = [x[:] for x in [[0] * totCol] * totRow]
    for i in range(len(ra)):
        row = (int)((dec[i] + 90)//rectH)
        col = (int)(ra[i]//rectW)
        rectColor[row][col] += Color[i]
        rectColorPointCount[row][col] += 1
    for row in range(totRow):
        for col in range(totCol):
            if rectColorPointCount[row][col] < whiteOutLimit:
                rectColor[row][col] = -1
            else:
                rectColor[row][col] = rectColor[row][col] / rectColorPointCount[row][col]
    np.savetxt("rectColor.cvs", rectColor, delimiter=",", fmt='%3.3f')
    if False:
        for row in range(totRow):
             for col in range(totCol):
                rectColorPointCount[row][col] = round((rectColorPointCount[row][col]), 2)
    np.savetxt("rectColorPointCount.cvs", rectColorPointCount, delimiter=",", fmt='%3.3i')
    for row in range(totRow):
        for col in range(totCol):
            wRA = col * rectW
            if wRA >= 180:
                wRA = wRA - 360
            W = np.array([np.radians(wRA), np.radians(wRA + rectW) + fillGapW])
            if rectColor[row][col] == -1:
                fillColor = 'purple'
            else:
                index = (int)((rectColor[row][col] - downFilter) * 255/ (upFilter-downFilter))
                fillColor = colorArray[index]
                ax.fill_between(W, np.radians(row * rectH - 90), np.radians((row + 1) * rectH - 90 + fillGapH), facecolor = fillColor)
                #print(index)

    #ax.fill_between(np.arange(np.radians(-170), np.radians(-150), 0.01), np.radians(-20), np.radians(20), facecolor = colorArray[255])


nameIND = 0
dateIND = 1
app_magIND = 2
abs_magIND = 3
raIND = 4
decIND = 5
zIND = 6
SNeTypeIND = 7
ebvIND = 8

totColNum = 9

dustA = 0 #dust value to adjust distance calculation

name, date, app_mag, abs_mag, ra, dec, z, SNeType, ebv = readIn(fileName)

ra = purifyValues(ra, True, False)
dec = purifyValues(dec, True, True)
raDegree = raToDegree(ra)
zeroDecUp = 2
zeroDecDown = -2
z = purifyValues(z, False, False)
app_mag = str2float(app_mag)
abs_mag = str2float(abs_mag)
ebv = str2float(ebv)
raRad, decRad = degreeToRadian(ra, dec)

maxZ, minZ = findMaxMin(z)
print("Total number of name data: " + str(len(name)))
print("number of data points between " + str(zeroDecDown) + " and " + str(zeroDecUp) + " dec: " + str(findZeroDec(dec)))
print("z value range from " + str(maxZ) + " to " + str(minZ))
ave_z, dist, v = plotLinearH()
print(ave_z)
Color = calculateColor(ave_z, dist, v)
figure=plot_mwd(ra, dec, Color, True)




for i in range(len(ra)):
    if ra[i]>180:
        ra[i] = ra[i]-360
#print(decRad)

local_H0=[]
for i in range((len(dist))):
    local_H0.append(v[i]/dist[i])
data_to_sample=[]






data_to_sample.append(np.radians(ra))
data_to_sample.append(np.radians(dec))


data_to_sample.append(local_H0)


def sampling(data,n_rows,n_cols):
    # data
    x = data[0]

    y = data[1]

    z = data[2]

    # prep
    x_i = 2*np.pi/n_cols # x increment
    y_i = np.pi/n_rows # y increment
    x_low = -np.pi# initial lower bound
    y_low = -np.pi/2 # initial lower bound

    # loop
    sampled_data = []
    for i in range(n_rows):
        x_low = -np.pi # initial lower bound
        for j in range(n_cols):
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
            x_low += x_i
            print('\rSampling Data: {:.1f}% done '.format((n_cols*i+j)*100/(n_rows*n_cols)),end='')
        y_low += y_i
    print('')

    # return data: [[rad,rad,rad,rad,val], . . . ]
    return sampled_data

data_sampled=sampling(data_to_sample,20,20)

def plot_rect(data):
    data = np.array(data)

    data = data.transpose()

    print(data)
    # TODO: still need to work on this
    # plot test: plots in mollweide projection with log scaling
    # increase niter for faster plot (plots every 'niter'th point
    plt.figure(figsize=(10,6))
    ax = plt.subplot(111,projection='mollweide')
    ax.scatter((data[0]+data[1])/2,(data[2]+data[3])/2,c=np.log10(data[4]),s=100)
    plt.draw()

    return 0
plot_rect(data_sampled)
plt.show()


#image= np.array(plt.imread('Figure_1.png',2))
#myplot=plt.imshow(image)
#plt.colorbar(myplot)
#plt.show()
