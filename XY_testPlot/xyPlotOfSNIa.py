import numpy as np
import matplotlib.pyplot as plt 
import csv

##
#@source: 1.Open Supernova Catalog
#   2.https://www.cfa.harvard.edu/supernova/README_archive
#@author: Xuyang Yu
#plot the distribution of Type Ia supernova on a XY-grid,
#where x-axis is the right ascension and y-axis is the declination
#domain: x is [0, 360], y is [-90, 90]

#Ia-norm    : spectroscopically/photometrically "normal" type-Ia supernovae
# Ia-91t     : type-Ia supernovae whose lightcurves/spectra resemble those of SN 1991T
# Ia-91bg    : type-Ia supernovae whose lightcurves/spectra resemble those of SN 1991bg
# Ia-pec     : all other type-Ia supernovae (e.g., 2000cx)
#

def hasColumn(row, index):
    if row[index] == '':
        return False
    return True

def readIn(fileName):
    name = []
    ra = []
    dec = []
    z = []
    SNeType = []
    csv_file = open(fileName, 'r')
    data_reader = csv.reader(csv_file, delimiter = ',')
    next(data_reader)
    for row in data_reader:
        if hasColumn(row, raIND) and hasColumn(row, decIND):
            name.append(row[nameIND])
            ra.append(row[raIND])
            dec.append(row[decIND])
            z.append(row[zIND])
            SNeType.append(row[SNeTypeIND])
    return name, ra, dec, z, SNeType

def str2float(data_list):
    for i in range(len(data_list)):
        data_list[i] = float(data_list[i])
    return data_list

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

fileName = 'OSCwithRA_DEC_z.csv'

nameIND = 0
raIND = 1
decIND = 2
zIND = 3
SNeTypeIND = 4

name, ra, dec, z, SNeType = readIn(fileName)

ra = purifyValues(ra, True, False)
dec = purifyValues(dec, True, True)
raDegree = raToDegree(ra)
zeroDecUp = 2
zeroDecDown = -2
z = purifyValues(z, False, False)

maxZ, minZ = findMaxMin(z)
print("Total number of data: " + str(len(name)))
print("number of data points between " + str(zeroDecDown) + " and " + str(zeroDecUp) + " dec: " + str(findZeroDec(dec)))
print("z value range from " + str(maxZ) + " to " + str(minZ))

plt.figure(figsize = (15,5))
plt.scatter(raDegree, dec, s = 1, c = 'blue')
plt.title('SNIa angular distribution from the OSC Data )')
plt.xlabel('Rigth ascension [degrees]')
plt.ylabel('declination [degrees]')

drawRectangles = False

if drawRectangles:
    for i in range(37):
        plt.axvline(i * 10, ymin = 0.02, ymax = 0.96, color='k')
    for i in range(19):
        plt.axhline(i * 10 - 90, xmin = 0.03, xmax = 0.97, color='k')

plt.show()

