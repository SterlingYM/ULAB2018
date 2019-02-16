import numpy as np
import matplotlib.pyplot as plt 
import csv

##
#@source: Hubble_Plot.py author: UK
#Xuyang's test python code
#add column host_name, EBV, and x_ray
#take average of all z_values (need to change)
#did not delete from duplicates
#change _err to _rej to better signify "rejected data" rather than "measurement error"
#

def isData(row):
    for i in range(len(row) - 1):
        if row[i] == '':
            return False
    return True

def readIn(fileName):
    name = []
    date = []
    app_mag = []
    abs_mag = []
    host_name = []
    z = []
    SNeType = []
    EBV = []
    err_row = []

    csv_file = open(fileName, 'r')
    data_reader = csv.reader(csv_file, delimiter = ',')
    next(data_reader)
    for row in data_reader:
        if isData(row):
            name.append(row[name_index])
            date.append(row[date_index])
            app_mag.append(row[app_mag_index])
            abs_mag.append(row[abs_mag_index])
            host_name.append(row[host_name_index])
            z.append(row[z_index])
            SNeType.append(row[type_index])   
            EBV.append(row[EBV_index])
        else:
            err_row.append(row)
    return name, date, app_mag, abs_mag, host_name, z, SNeType, EBV, err_row         

def purify_z(Z):
    for i in range(len(Z)):
        strList = Z[i].split(',')
        z_a = 0.0
        for j in range(len(strList)):
            z_a = z_a + (float)(strList[j])
        Z[i] = (str)(z_a / len(strList))
    return Z

def str2float(data_list):
    for i in range(len(data_list)):
        data_list[i] = float(data_list[i])
    return data_list

def dist_velocity(name, app_mag, abs_mag, z):
    dist = []
    dist_rej = []
    v = []
    v_rej = []
    approved_SNe_data = []
    rejected_SNe_data = []

    for i in range(len(z)):
        x = 10 ** ((app_mag[i] - abs_mag[i] + 5) / 5)
        y = z[i] * c
        if y not in v:
            dist.append(x)
            v.append(y)
            approved_SNe_data.append([name[i], app_mag[i], abs_mag[i], z[i]])
        else:
            rejected_SNe_data.append([name[i],app_mag[i],abs_mag[i],z[i]])
            dist_rej.append(x)
            v_rej.append(y)
    return approved_SNe_data, dist, v, rejected_SNe_data, dist_rej, v_rej

c = 3.0 * 10 ** 5 #unit [km/s]

filename = 'TestData1.csv'

name_index = 0
date_index = 1
app_mag_index = 2
abs_mag_index = 3
host_name_index = 4
z_index = 5
type_index = 6
EBV_index = 7

name, _, app_mag, abs_mag, _, z, _, _, err_row = readIn(filename)

app_mag = str2float(app_mag)
abs_mag = str2float(abs_mag)
z = purify_z(z)
z = str2float(z)  

approved_data, dist, v, rejected_data, dist_rej, v_rej = dist_velocity(name, app_mag, abs_mag, z)

plt.figure(figsize = (15,5))
plt.scatter(dist, v, s = 5, c = 'blue')
plt.show()

print("Total number of data: " + str(len(v)))
