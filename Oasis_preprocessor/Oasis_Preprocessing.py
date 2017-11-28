
# coding: utf-8

# """
# Last Revised 11/20/2017
# Jordan Wilson
# USGS
# Missouri Water Science Center
# jlwilson@usgs.gov
# 
# Combines multiple resistivity and water-quality profiles in the form of semicolon and comma delimited .txt files into .csv files for import into Oasis. 
# 
# Tkinter code was modified from DSWallace integrateGeophysicsQW.py script. 
# """
#%%
from Tkinter import *
import pandas as pd
from Tkinter import Tk
from tkFileDialog import asksaveasfilename
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename
import tkMessageBox
import os
from shapely.geometry import Point
import geopandas as gp
import numpy as np
import shutil, glob
import fiona

#%%
# Defining bandpass filter to filter resistivity data channels
def band_pass(df,column, min, max):
    df[column+'_bandpass']=df[column]
    df[column+'_bandpass'][df[column]<min] = np.nan
    df[column+'_bandpass'][df[column]>max] = np.nan

#%%
# Defining depth filter
def depth_filt(df, column, offset, factor):
    df[column+'_filt']=df[column]
    df[column+'_filt'][df[column]<offset+factor]=np.nan
    
#%%
# Defining rolling average filter
def rolling_avg(df, column1, column2, width):
    df[column1+'_rollavg']=df[column2]
    df[column1+'_rollavg']= pd.rolling_mean(df[column1+'_rollavg'],window=width)

    
#%%
Tk().withdraw() 
tkMessageBox.showinfo("Directions", "For this script to work, you must have all resistivity .txt files that you want to combine in one folder. All QW .csv files must also be in one folder. It may be the same folder.")

# Hard-coded for now, tired of selecting.  Will remove when finalized
res_folder = "D:\\Mississippi Alluvial Plain\\2018\\Waterborne Resistivity Scripts\\TestFiles\\FromWilson\\Pre_Oasis\\Resistivity"
wq_folder = "D:\\Mississippi Alluvial Plain\\2018\\Waterborne Resistivity Scripts\\TestFiles\\FromWilson\\Pre_Oasis\\QW"
ini_file = "D:\\Mississippi Alluvial Plain\\2018\\Waterborne Resistivity Scripts\\TestFiles\\FromWilson\\Pre_Oasis\\Resistivity\\Floodway_Inflatable_10m_extension_06142017.ini"
directory = "D:\\Mississippi Alluvial Plain\\2018\\Waterborne Resistivity Scripts\\TestFiles\\FromWilson\\Pre_Oasis\\TestOutput"

"""
# Resistivity files
res_folder = askdirectory(title="Select folder that contains all raw resistivity files for processing...")  # show an "Open" dialog box and return the path to the selected file
if not glob.glob('{}/*.txt'.format(res_folder)):
    tkMessageBox.showerror("FILE ERROR", "No resistivity files contained within folder or incorrect format")
    exit()

# Water Quality Files
wq_folder = askdirectory(title="Select folder that contains all raw water-quality data for processing...")
if not glob.glob('{}/*.csv'.format(wq_folder)):
    tkMessageBox.showerror("FILE ERROR", "No water quality files contained within folder or incorrect format")
    exit()

# Initialization File
ini_file = askopenfilename(title="Select ini file used to collect the resistivity data",
                           filetypes=[("INI Files", "*.ini")])
if not ini_file:
    tkMessageBox.showerror("FILE ERROR", "No INI file selected")
    exit()

# Save File Location
directory = askdirectory(title="Select directory to save the reordered resistivity and water-quality data")
"""
"""
# %% NHD geodatabase Location
gdb_dir = askdirectory(title="Select directory where NHD geodatabase is located")
try:
    streams = gp.read_file(gdb_dir)
except:
    tkMessageBox.showerror("FILE ERROR", "No NHD geodatabase selected")
    exit()
stream_df=pd.DataFrame(streams)
#%%
str(stream_df.ix[3,'geometry'])
#%%
for x in range(0,len(stream_df.FCode)):
    stream_df.ix[x,'geo']=str(stream_df.ix[x,'geometry'])
#%%

#Rename files from upstream to downstream with indicator of direction '_up' and '_down'

#Ideas: tkinter prompt to select if files are upstream or downstream 
#%%

"""
#%%

# Create crosstab table and export as a .txt file showing old and new names

#%%
path = directory + r'/Raw_Data_Renamed'

try:
    os.makedirs(os.path.join(path))
    print('Raw data folder created')
except:
    pass  # This makes me nervous- what errors are you avoiding? dsw 20171128

#%%

# Export renamed files to 'Raw_Data_Renamed' directory
    
#%%
# Preprocessing Resistivity Data
print('Aggregating raw data files')
outfilename="{}/all.txt".format(res_folder)

# Copy resistivity data into a single file
with open(outfilename, 'wb') as outfile:
    for filename in glob.glob('{}/*.txt'.format(res_folder)):
        if filename == outfilename:
            continue
        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)

# Reformat "all" delimiters to semicolons
replacements = {',':';'}
lines = []
with open('{}/all.txt'.format(res_folder)) as infile:
    for line in infile:
        for src, target in replacements.iteritems():
            line = line.replace(src, target)
        lines.append(line)
with open('{}/all.txt'.format(res_folder), 'w') as outfile:
    for line in lines:
        outfile.write(line)

print('Importing resistivity data')
importfile = pd.read_csv('{}/all.txt'.format(res_folder), sep=';', index_col=False, skiprows=1,
                         names=["Distance", "Depth", "Rho 1", "Rho 2", "Rho 3", "Rho 4", "Rho 5", "Rho 6", "Rho 7",
                                "Rho 8", "Rho 9", "Rho 10", "C1", "C2", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8",
                                "P9", "P10", "P11", "Latitude", "Longitude", "In_p", "In_n", "V1_p", "V1_n", "V2_p",
                                "V2_n", "V3_p", "V3_n", "V4_p", "V4_n", "V5_p", "V5_n", "V6_p", "V6_n", "V7_p", "V7_n",
                                "V8_p", "V8_n", "V9_p", "V9_n", "V10_p", "V10_n", "GPSString", "UTC", "Latitude2", "D1",
                                "Longitude2", "D2", "Fix Quality", "Satellites", "HDOP", "Altitude", "D3",
                                "Height of Geoid"],
                         usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                  24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
                                  45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60])

print('Processing resistivity data')
data_cols = ('Ohm_m',
             'Cor_Dist',
             'Cor_Depth',
             'Final_Rho_1',
             'Final_Rho_2',
             'Final_Rho_3',
             'Final_Rho_4',
             'Final_Rho_5',
             'Final_Rho_6',
             'Final_Rho_7',
             'Final_Rho_8',
             'Final_Rho_9',
             'Final_Rho_10',
             'Lat',
             'Lon',
             'Final_Altitude',
             'File')
# Add additional data columns listed above and remove unwanted data columns
for x in data_cols:
    importfile[x]=np.nan
importfile.drop(['D1','D2','D3', 'File'], inplace=True, axis=1)
importfile.insert(0,'File',"")

# Demarcate which file data came from and remove headers from combined data --> NOTE: THIS TAKES A WHILE TO RUN
j=1
for x in range(0,len(importfile.Distance)):
    if importfile.ix[x,"Distance"]=="Distance":
        j+=1
    elif importfile.ix[x,"Distance"]!="Distance":
        importfile.ix[x,"File"]=j
    else:
        break
# Removes bad data files
importfile = importfile.loc[importfile.GPSString == 'GPGGA', :]

"""
# Propose dropping this section for simplification - can't see that it does anything
importfile = importfile.loc[(importfile.GPSString=='GPGGA')|(importfile.Distance=='Distance'),:]


importfile['Counter']=importfile.index.values+1
counter = importfile['Counter']
importfile.drop(['Counter'], axis=1, inplace=True)
importfile.insert(0,'Counter',counter)

importfile = importfile.loc[importfile.Distance.str.contains("Distance")==False,:]
"""
importfile.set_index([range(0,len(importfile.Distance))], inplace=True)

# Reformat Latitude and Longitude to decimal degrees
importfile['Lat1'] = importfile['Latitude'].str[0:2]
importfile['Lat2'] = importfile['Latitude'].str[2:]
importfile['Lat2'] = importfile['Lat2'].astype(float)
importfile['Lat1'] = importfile['Lat1'].astype(float)
importfile['Lat'] = importfile.Lat1+importfile.Lat2/60

importfile['Lon1'] = importfile['Longitude'].str[0:3]
importfile['Lon2'] = importfile['Longitude'].str[3:]
importfile['Lon2'] = importfile['Lon2'].astype(float)
importfile['Lon1'] = importfile['Lon1'].astype(float)
importfile['Lon'] = importfile.Lon1-importfile.Lon2/60

importfile.drop(['Lat1','Lat2','Lon1','Lon2'], axis=1, inplace=True)
# importfile.drop('Counter', axis=1, inplace=True)  # Propose dropping this code

# Reformat numbers in the file to float
for col in importfile.columns[1:]: 
    try:
        importfile[col] = importfile[col].astype(float)
    except:
        pass

# Creates new column with new incremental file count
importfile['File1']="1"
j=1
for x in range(1,len(importfile.Distance)):
    if importfile.ix[x,"File"]==importfile.ix[x-1,"File"]:
        importfile.ix[x,"File1"]=importfile.ix[x-1,"File1"]
        
    elif importfile.ix[x,"File"]!=importfile.ix[x-1,"File"]:
        j+=1
        importfile.ix[x,"File1"]=j
    else:
        break
file1 = importfile['File1']
importfile.drop(['File1', 'File'], axis=1, inplace=True)
importfile.insert(0, 'File', file1)

#%%
# Calculating the distance from UTM coordinates
importfile['Cum_dist']=0

for i in range(1,len(importfile['Distance'])):
    importfile.ix[i,'Cor_Dist']=np.sqrt(np.square(importfile.ix[i,'X_UTM']-importfile.ix[i-1,'X_UTM'])+np.square(importfile.ix[i,'Y_UTM']-importfile.ix[i-1,'Y_UTM']))
    importfile.ix[i,'Cum_dist']=importfile.ix[i-1,'Cum_dist']+importfile.ix[i,'Cor_Dist']

# Here is where we might search for gaps in the data...

#%%
# Import INI file
        
ini = pd.read_csv(ini_file, index_col=None,sep='=')
depthoffset = float(ini.ix['DepthOffset','[SwitchPro]'])

#%%
# Applying the bandpass filter and rolling average
for x in range(1,11):
    band_pass(importfile,'Rho {}'.format(x),0,250)
    rolling_avg(importfile, 'Rho {}'.format(x), 'Rho {}_bandpass'.format(x), 20)

#%%
# Applying the depth filter
depth_filt(importfile, 'Depth', depthoffset, 0.01)
rolling_avg(importfile, 'Depth', 'Depth_filt'.format(x), 20)

#%%
# Filtering Altitude via rolling median filter
importfile['Altitude_filt'] = pd.rolling_median(importfile['Altitude'],window=10)
rolling_avg(importfile, 'Altitude', 'Altitude_filt', 20)
importfile['Altitude_rollavg']=importfile['Altitude_rollavg'].round(1)

#%%
# Converting WGS 84 coordinates to UTM 15N coordiantes
geometry = [Point(xy) for xy in zip(importfile.Lon, importfile.Lat)]
crs=None
importfile = gp.GeoDataFrame(importfile, crs=crs, geometry=geometry)
importfile.crs = {'init' :'epsg:4326'}
importfile = importfile.to_crs({'init': 'epsg:32615'})

def getXY(pt):
    return (pt.x, pt.y)
centroidseries = importfile['geometry'].centroid
x,y = [list(t) for t in zip(*map(getXY, centroidseries))]

importfile['X_UTM']=x
importfile['Y_UTM']=y

#%%
#Replacing all NaNs with "*" and dropping geometry column
importfile.fillna('*', inplace=True)
importfile.drop('geometry', axis=1, inplace=True)

#%%

saveRes = asksaveasfilename(defaultextension='.csv',title="Designate resitivity csv name and location", filetypes=[('csv file', '*.csv')])
importfile.to_csv(saveRes, index=False)
print('Resistivity data exported')




#%%
# Preprocessing QW Data

from os import listdir
from os.path import isfile, join

mypath=wq_folder

qwfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

j=qwfiles[0]

print('Importing water quality data')
data = pd.read_csv('{}/{}'.format(wq_folder,j),sep=',',skiprows=12,index_col=False,engine='python',encoding='utf-16',
                   names=["Date","Time","°C","mmHg","DO %","SPC-uS/cm","C-uS/cm","ohm-cm","pH","NH4-N mg/L",
                          "NO3-N mg/L","Cl mg/L","FNU","TSS mg/L","DEP m","ALT m","Lat","Lon"])

qwdata=data
qwdata['File']=1
z=2
try:
    for x in qwfiles[1:]:
        data = pd.read_csv('{}/{}'.format(wq_folder,x),sep=',',skiprows=12,index_col=False,engine='python',
                           encoding='utf-16', names=["Date","Time","°C","mmHg","DO %","SPC-uS/cm","C-uS/cm",
                                                     "ohm-cm","pH","NH4-N mg/L","NO3-N mg/L","Cl mg/L","FNU",
                                                     "TSS mg/L","DEP m","ALT m","Lat","Lon"])
        qwdata['File']=z
        z+=1
        qwdata=qwdata.append(data)

    print('Water quality data imported')
except:
    pass
print('Processing water quality data')
qwdata.dropna(axis=1, how='all', inplace=True)

qwdata.rename(columns={'°C':'Temp_C', 'SPC-uS/cm': 'SPC_mscm','C-uS/cm':'Cond_mscm','ohm-cm':'Res_ocm','ALT m':'Alt_m'}, inplace=True)
try:
    qwdata = qwdata.loc[qwdata.Res_ocm!="    +++++",:]
except:
    pass
qwdata['Res_ocm'] = qwdata['Res_ocm'].astype('float')
qwdata['Ohm_m']=qwdata['Res_ocm']/100

#%%
# Applying a rolling average on resistivity
rolling_avg(qwdata, 'Ohm_m', 'Ohm_m', 20)

#%%
for col in qwdata.columns[1:]:                  
    try:
        qwdata[col] = qwdata[col].astype(float)
    except:
        pass

#%%
# Converting WGS 84 coordinates to UTM 15N coordiantes
geometry = [Point(xy) for xy in zip(qwdata.Lon, qwdata.Lat)]
crs=None
qwdata = gp.GeoDataFrame(qwdata, crs=crs, geometry=geometry)
qwdata.crs = {'init' :'epsg:4326'}
qwdata = qwdata.to_crs({'init': 'epsg:32615'})

def getXY(pt):
    return (pt.x, pt.y)
centroidseries = qwdata['geometry'].centroid
x,y = [list(t) for t in zip(*map(getXY, centroidseries))]

qwdata['X_UTM']=x
qwdata['Y_UTM']=y

#%%
# Filling all NaNs with '*'
qwdata.drop('geometry',axis=1,inplace=True)
qwdata.fillna('*', inplace=True)

#%%
saveQW = asksaveasfilename(defaultextension='.csv',title="Designate water quality csv name and location", filetypes=[('csv file', '*.csv')])
qwdata.to_csv(saveQW, index=False)
print('Water quality data exported!')
