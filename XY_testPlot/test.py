
import numpy as np
import matplotlib.pyplot as plt 
import csv

def isData(row):
	for i in range(len(row)):
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
    x_ray = []
    err_row = []

    csv_file = open(fileName, 'r')
    data_read = csv.reader(csv_file, delimiter = ',')
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
            x_ray.append(row[x_ray_index]) 		

