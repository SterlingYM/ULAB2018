#!/usr/bin/env python
# coding: utf-8

# In[41]:


import numpy as np
import matplotlib.pyplot as plt
import csv

####### functions #######

def isData2use(row,limitType):
    # Returns True if data is to be used
    for i in range(len(row)):
        #return false if any type of data is missing
        if row[i] == '':
            return False
    # if none of data is empty, data is goot to go
    # check 'Type' if limitType == True:
    if limitType:
        return row[type_index]=='Ia'
    return True

def readIn(filename,limitType):
    # This function reads in csv file and store apparent/absolute magnitudes and z-value.
    # If there is a data with missing values, that datapoint(=row) is stored in err_row
    # Note that these returned values are still string
    name = []
    date = []
    app_mag = []
    abs_mag = []
    z = []
    SNeType = []
    err_row = []
    
    csv_file = open(filename,'r') #readonly
    data_reader = csv.reader(csv_file, delimiter = ',')
    next(data_reader) #skip header
    for row in data_reader:
        if isData2use(row,limitType):
            name.append(row[name_index])
            date.append(row[date_index])
            app_mag.append(row[app_mag_index])
            abs_mag.append(row[abs_mag_index])
            z.append(row[z_index])
            SNeType.append(row[type_index])
        else:
            err_row.append(row)
            
    return name, date, app_mag, abs_mag, z, SNeType, err_row

def purify_z(Z):
    # This data sometimes contain two values for z separated by ','.
    # Those values are usually equal. 
    # This function 'purifies' those duplicated values into single value.
    # Returned values are still string, so it must be converted to float afterward.
    for i in range(len(Z)):
        str = Z[i]
        if str.find(','):
            Z[i] = str[:str.find(',')]
    return Z

    
def str2float(data_list):
    # This function takes any string list and convert into float list.
    for i in range(len(data_list)):
        data_list[i] = float(data_list[i])
    return data_list


def dist_velocity(name,app_mag,abs_mag,z):
    # This calculated data contain lots of groups of points that have same velocity values.
    # Having exactly same velocity is unnatural (probably due to low precision).
    # This function takes apparent/absolute magnitudes and redshift value
    # and calculates distance and velocity.
    # If there is already the same value in list, it will be stored in _err lists.
    dist = []
    dist_err = []
    v = []
    v_err = []
    approved_SNe_data = []
    rejected_SNe_data = []
    # data list formatted in [name,app_mag,abs_mag,z]
    
    for i in range(len(z)):
        x = 10**((app_mag[i] - abs_mag[i] + 5)/5) # distance modulus
        y = z[i]*c     # radial velocity
        if y not in v: # reject if y-data matches with existing
            dist.append(x) 
            v.append(y)
            approved_SNe_data.append([name[i],app_mag[i],abs_mag[i],z[i]])
        else:
            # find existing duplicate, add to error list, and delete from approved list
            dup_index = v.index(y)
            rejected_SNe_data.append(approved_SNe_data[dup_index])
            dist_err.append(dist[dup_index])
            v_err.append(v[dup_index])
            del approved_SNe_data[dup_index]
            del dist[dup_index]
            del v[dup_index]
            
            # add new duplicate to error list
            rejected_SNe_data.append([name[i],app_mag[i],abs_mag[i],z[i]])
            dist_err.append(x)
            v_err.append(y)
            
    return approved_SNe_data, dist, v, rejected_SNe_data, dist_err, v_err




####### constants #######

# SOL: unit of [km/s] is typically used
c = 3.0 * 10 ** 5

# csv file: This file should contain: [Name, Date, m, M, z, Type]
filename = 'TheOpenSupernovaCatalog.csv'

#index of each data type in 'row'. Change this if data format is different
name_index = 0
date_index = 1
app_mag_index = 2
abs_mag_index = 3
z_index = 4
type_index = 5

# change this to 'False' for all-type plot (not Ia-only)
limitType = True




####### code body #######

# (1) read csv file
# not importing discovered date and SNe type
# changing '_' to proper variable names will allow importing those values
name, _, app_mag, abs_mag, z, _, err_row = readIn(filename,limitType)

# (2) change data values into usable form
app_mag = str2float(app_mag)
abs_mag = str2float(abs_mag)
z = purify_z(z)
z = str2float(z)    

# (3) calculate values to plot
approved_data, dist, v, rejected_data, dist_err, v_err = dist_velocity(name, app_mag, abs_mag, z)

# (4) plot
plt.figure(figsize = (15,5))
plt.scatter(dist_err,v_err,s=5,c='orange')
plt.scatter(dist,v,s=5,c='blue')
plt.title('The Open Supernova Catalog Data (Type Ia)')
plt.xlabel('distance [pc]')
plt.ylabel('velocity [km/s]')
plt.show()

plt.figure(figsize = (15,5))
plt.scatter(dist,v,s=5,c='blue')
plt.title('The Open Supernova Catalog Data (Type Ia, picked)')
plt.xlabel('distance [pc]')
plt.ylabel('velocity [km/s]')
plt.show()

# (5) description
print('Total number of data: {}'.format(len(z)+len(err_row)))
print('Number of data with missing values: {}'.format(len(err_row)))
print('Number of data with exact duplicate of v value (orange plot): {}'.format(len(v_err)))
print('Number of \'good\' data (blue plot): {}'.format(len(v)))
print('This descrete values of velocity is considered due to low sig-fig of z-data.')

print('\nSample of \'high precision\' data (blue-plot):')
print('{:15}{:10}{:10}{:11}'.format('Name','m','M','z'))
for i in range(10): print('{:11}{:10.3f}{:10.3f}\t   {:<14}'.format(*approved_data[i])) 
    
print('\nSample of \'low precision\' data (orange-plot):')
print('{:15}{:10}{:10}{:11}'.format('Name','m','M','z'))
for i in range(10): print('{:11}{:10.3f}{:10.3f}\t   {:<14}'.format(*rejected_data[i])) 
    
print('\nFrom this data analysis, it is clear that precision of spectroscopic data is critical.')


# Part 2: fitting model with data
# 
# The word "fitting" is often used to draw a line on the data.
# this "fitting" requires two steps:
#
#   (1) construct a model
#       in cosmology, 'model' refers to a set of theoretical equations.
#       'model' is coded as a function, and this model uses a very early model
#       with assumption of constant H.
#
#   (2) find best-fit parameters
#       once model is prepared, we use a fitting function that automatically
#       finds the best-fit value for parameters.
#       see documentation of scipy.optimize.
#       astropy has also fitter, and it is often used for more complicated model.


import scipy.optimize as fitter

def my_model(x,H):
    return H*x

def fitting(model,xdata,ydata):
    # data fitting
    par0 = np.array([1])
    par, cov = fitter.curve_fit(model, xdata, ydata, par0, absolute_sigma=True)
    
    #test
    #print(par[0])

    # plot
    x = np.linspace(0,max(xdata),10**4)
    y_fitted = my_model(x,par[0])
    
    plt.figure(figsize = (15,5))
    plt.plot(x,y_fitted,c='orange')
    plt.scatter(xdata,ydata,s=5)
    plt.title("Fitted function")
    plt.legend(['da/dt = {:2f}a'.format(par[0])])
    plt.xlabel('distance [pc]')
    plt.ylabel('velocity [km/s]')
    plt.show()

    # result
    print("Hubble Constant H0 = {}".format(par[0]))


print("\n\nPart 2: fitting")
fitting(my_model,dist,v)
