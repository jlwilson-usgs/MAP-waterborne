
# coding: utf-8

# In[31]:

"""
# Last Revised 1/18/2018
#
# Jordan Wilson
# USGS
# Missouri Water Science Center
# jlwilson@usgs.gov
#
# David Wallace
# USGS
# Texas Water Science Center
# dswallace@usgs.gov
#
# Combines multiple resistivity and water-quality profiles in the form of semicolon and comma delimited .txt files and exports .csv files for import into Oasis and .shp for GIS.
#
#
"""
#%%
import easygui as eg
import sys
import pandas as pd
from Tkinter import *
from Tkinter import Tk
from tkFileDialog import asksaveasfilename
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename
import tkMessageBox
import tkSimpleDialog
import os
from shapely.geometry import Point
import geopandas as gp
import numpy as np
import glob
from math import radians, cos, sin, asin, sqrt
from shutil import copyfile
import logging
import traceback
import datetime
import warnings
import scipy.stats

#%%
def workbench():
    global out
    eg.msgbox("For this script to work, the final Rho columns must be named 'Final_Rho_1,Final_Rho_2,...,Final_Rho_n', Latitude and Longitude columns must be named 'Lat' and 'Lon', and the corrected distance and depth columns should be named 'Cor_Dist' and 'Cor_Depth', and the QW columns containing resisitivty values from the QW meter should be named 'Ohm_m'")
    infile = eg.fileopenbox(title="Select Oasis output csv file for processing")

    try:
        print(userRiverName)
    except:
        userRiverName = tkSimpleDialog.askstring("River Reach", "Please enter the name of the river reach...", initialvalue="RIVER")

    try:
        data = pd.read_csv(infile)
        logging.info("Oasis file read\n")
    except:
        logging.info("Oasis file not read\n")

    try:
        data['File_w']=data['Line'].str.split('L',1)
        data['File']=''
        data['File']=data.apply(lambda row: row['File_w'][1],axis=1)
        data['File']=data['File'].astype('int')

    except:

        try:
            data['File']=data['File'].astype('int')

        except KeyError,e:
            logging.info('The following column is missing in the input file: %s. Check to make sure all required columns are present.\n' % str(e))

    out = pd.DataFrame()

    data.rename(columns={'File':'Profile'}, inplace=True)

    headers=('/Water_Res','Cor_Dist','Cor_Depth','Rho_1','Rho_2','Rho_3','Rho_4','Rho_5','Rho_6','Rho_7','Rho_8','Rho_9','Rho_10','C1','C2','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','Lat','Lon','Final_Altitude','Profile')
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
                 'C1','C2','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11',
                 'Lat',
                 'Lon',
                 'Final_Altitude',
                 'Profile')

    zipped = zip(headers, data_cols)
    try:
        for x,y in reversed(zipped):
            out.insert(0,x,data['{}'.format(y)].values)
        logging.info("File formatted\n")
        outfile = eg.filesavebox(title="Save processed file as...",default='{}_WorkbenchImport.csv'.format(userRiverName),filetypes=['*.csv'])
        out.to_csv(outfile,index=False)
        logging.info("Output csv file saved\n")
        logging.info("Workbench Preprocessor Finished\n")
    except KeyError, e:
        logging.info("Missing column",'The following column is missing in the input file: %s. Check to make sure all required columns are present.\n' % str(e))
    #%%
def workbench_checks():
    for x in range(0,len(out.Profile)):

        # Line-to-Line Continuity Check
        if out.ix[x,"Profile"]==out.ix[x-1,"Profile"]:
            pass

        elif out.ix[x,"Profile"]!=out.ix[x-1,"Profile"]:
            for s in range(1,11):
                if np.abs((out.ix[x,"Rho_{}".format(s)]-out.ix[x-1,"Rho_{}".format(s)])/(out.ix[x,"Rho_{}".format(s)]+out.ix[x-1,"Rho_{}".format(s)])/2*100)>=50:
                    tkMessageBox.showwarning("WARNING", "Large relative percent difference (>50%) in Rho {} on Line {} / row {}".format(s,out.ix[x,"Profile"],x+2))
                    logging.warning("Line-to-Line Disconitnuity (>50%) in Rho {}".format(s))
                else:
                    pass
            if np.abs((out.ix[x,"/Water_Res".format]-out.ix[x-1,"/Water_Res"])/(out.ix[x,"/Water_Res"]+out.ix[x-1,"/Water_Res"])/2*100)>=50:
                    tkMessageBox.showwarning("WARNING", "Large relative percent difference (>50%) in Water_Res on Line {} / row {}".format(out.ix[x,"Profile"],x+2))
                    logging.warning("Line-to-Line Disconitnuity (>50%) in Water_Res")
            else:
                pass
        else:
            pass
    logging.info("Continuity check finished\n")

# Line-to-Line Min/Max Check
    mult=2 #Mulitiple of the inner quartile range (1st and 3rd) that defines outliers

    # Resistivity checks
    for x in out['Profile'].unique():
        for s in range(1,11):
            if np.max(out['Rho_{}'.format(s)][out['Profile']==x])>mult*(scipy.stats.mstats.mquantiles(out['Rho_{}'.format(s)][out['Profile']==x])[2]):
                tkMessageBox.showwarning("WARNING", "Maximum Rho_{} in row {} larger than {} times the 3rd quantile of Line {}".format(s,out.index[out['Rho_{}'.format(s)]==np.max(out['Rho_{}'.format(s)][out['Profile']==x])][0]+2,mult,x))
                logging.warning('Maximum Rho_{} in row {} larger than {} times the 3rd quartile of Line {}'.format(s,out.index[out['Rho_{}'.format(s)]==np.max(out['Rho_{}'.format(s)][out['Profile']==x])][0]+2,mult,x))
            elif np.min(out['Rho_{}'.format(s)][out['Profile']==x])<1/mult*(scipy.stats.mstats.mquantiles(out['Rho_{}'.format(s)][out['Profile']==x])[0]):
                tkMessageBox.showwarning("WARNING", "Minimum Rho_{} in row {} smaller than 1/{} times the 1st quartile of Line {}".format(s,out.index[out['Rho_{}'.format(s)]==np.min(out['Rho_{}'.format(s)][out['Profile']==x])][0]+2,mult,x))
                logging.warning("Minimum Rho_{} in row {} smaller than 1/{} times the 1st quartile of Line {}".format(s,out.index[out['Rho_{}'.format(s)]==np.min(out['Rho_{}'.format(s)][out['Profile']==x])][0]+2,mult,x))
            else:
                pass


    # QW checks
    for x in out['Profile'].unique():
        if np.max(out['/Water_Res'][out['Profile']==x])>mult*(scipy.stats.mstats.mquantiles(out['/Water_Res'][out['Profile']==x])[2]):
            tkMessageBox.showwarning("WARNING", "Maximum Water_Res in row {} larger than {} times the 3rd quantile of Line {}".format(out.index[out['/Water_Res']==np.max(out['/Water_Res'][out['Profile']==x])][0]+2,mult,x))
            logging.warning('Maximum Water_Res in row {} larger than {} times the 3rd quartile of Line {}'.format(out.index[out['/Water_Res']==np.max(out['/Water_Res'][out['Profile']==x])][0]+2,mult,x))
        elif np.min(out['/Water_Res'][out['Profile']==x])<1/mult*(scipy.stats.mstats.mquantiles(out['/Water_Res'][out['Profile']==x])[0]):
            tkMessageBox.showwarning("WARNING", "Minimum Water_Res in row {} smaller than 1/{} times the 1st quartile of Line {}".format(s,out.index[out['/Water_Res']==np.min(out['/Water_Res'][out['Profile']==x])][0]+2,mult,x))
            logging.warning("Minimum Water_Res in row {} smaller than 1/{} times the 1st quartile of Line {}".format(out.index[out['/Water_Res']==np.min(out['/Water_Res'][out['Profile']==x])][0]+2,mult,x))
        else:
            pass

    # Altitude checks
        if np.max(out['Final_Altitude'][out['Profile']==x])>mult*(scipy.stats.mstats.mquantiles(out['Final_Altitude'][out['Profile']==x])[2]):
            tkMessageBox.showwarning("WARNING", "Maximum Altitude in row {} larger than {} times the 3rd quantile of Line {}".format(out.index[out['Final_Altitude']==np.max(out['Final_Altitude'][out['Profile']==x])][0]+2,mult,x))
            logging.warning('Maximum Altitude in row {} larger than {} times the 3rd quartile of Line {}'.format(out.index[out['Final_Altitude']==np.max(out['Final_Altitude'][out['Profile']==x])][0]+2,mult,x))
        elif np.min(out['Final_Altitude'][out['Profile']==x])<1/mult*(scipy.stats.mstats.mquantiles(out['Final_Altitude'][out['Profile']==x])[0]):
            tkMessageBox.showwarning("WARNING", "Minimum Altitude in row {} smaller than 1/{} times the 1st quartile of Line {}".format(out.index[out['Final_Altitude']==np.min(out['Final_Altitude'][out['Profile']==x])][0]+2,mult,x))
            logging.warning("Minimum Altitude in row {} smaller than 1/{} times the 1st quartile of Line {}".format(out.index[out['Final_Altitude']==np.min(out['Final_Altitude'][out['Profile']==x])][0]+2,mult,x))
        else:
            pass

    logging.info("Min/max check finished\n")

# Blank Cells Check
    if out.isnull().values.any():
        nulls = out.isnull()
        for t in range(0,len(out.columns)):
            for x in range(0,len(out.Profile)):
                if nulls.ix[x,t]:
                    tkMessageBox.showwarning("WARNING", "Blank {} cell in row {}".format(out.columns.values[t],x+2))
                    logging.warning('Blank {} cell in row {}'.format(out.columns.values[t],x+2))
                else:
                    pass
    else:
        pass
    logging.info("Blank cell check finished\n")

#%%
def oasis():
    #%%
    global userRiverName, importfile, importfile1
    # Supressing depreciation warning from output

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
        df[column1+'_rollavg']= df[column1+'_rollavg'].rolling(width, min_periods=1).mean()

    #%%
    # Defining the rolling median filter
    def rolling_median(df, column1, column2, width):
        df[column1+'_rollmed']=df[column2]
        df[column1+'_rollmed']= df[column1+'_rollmed'].rolling(width, min_periods=1).median()

    def haversine(lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a))
        r = 6371  # Radius of earth in kilometers. Use 3956 for miles
        return c * r

    # %% -----------------------------------------------------------------------------------------------------------------
    Tk().withdraw()
    tkMessageBox.showinfo("Directions", "For this script to work, you must have all resistivity .txt files that you want to combine in one folder. All QW .csv files must also be in one folder. It may be the same folder.")

    userRiverName = tkSimpleDialog.askstring("River Reach", "Please enter the name of the river reach...",
                                             initialvalue="RIVER")

    # Resistivity files
    res_folder = askdirectory(title="Select folder that contains all raw resistivity files for processing...")  # show an "Open" dialog box and return the path to the selected file
    if not glob.glob('{}/*.txt'.format(res_folder)):
        tkMessageBox.showerror("FILE ERROR", "No resistivity files contained within folder or incorrect format")
        logging.error("No text files found within selected resistivity folder\n")
        exit()

    # Water Quality Files
    wq_folder = askdirectory(title="Select folder that contains all raw water-quality data for processing...", initialdir=res_folder)
    if not glob.glob('{}/*.csv'.format(wq_folder)):
        tkMessageBox.showerror("FILE ERROR", "No water quality files contained within folder or incorrect format")
        logging.error("No csv files found within the selected water quality folder\n")
        exit()

    # Initialization File
    ini_file = askopenfilename(title="Select ini file used to collect the resistivity data",filetypes=[("INI Files", "*.ini")], initialdir=res_folder)
    if not ini_file:
        tkMessageBox.showerror("FILE ERROR", "No INI file selected")
        logging.error("No INI file selected by the user\n")
        exit()

    # Save File Location
    directory = askdirectory(title="Select directory to save the reordered resistivity and water-quality data", initialdir=res_folder)

    # %% -----------------------------------------------------------------------------------------------------------------
    path = directory + r'/Raw_Data_Renamed'

    try:
        os.makedirs(os.path.join(path))
        print('Raw data folder created')
        logging.info("Raw data folder created\n")
    except:
        logging.info("Unable to create raw data folder or folder already exists\n")
        pass  # This makes me nervous- what errors are you avoiding? dsw 20171128

    # %% -----------------------------------------------------------------------------------------------------------------
    print("Verifying continuity of surveys")
    logging.info("Verifying continuity of surveys\n")
    # Correct names of columns in resistivity raw data files
    colNames = ["Distance", "Depth", "Rho 1", "Rho 2", "Rho 3", "Rho 4", "Rho 5", "Rho 6", "Rho 7", "Rho 8", "Rho 9",
                "Rho 10", "C1", "C2", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "Latitude",
                "Longitude", "In_p", "In_n", "V1_p", "V1_n", "V2_p", "V2_n", "V3_p", "V3_n", "V4_p", "V4_n", "V5_p", "V5_n",
                "V6_p", "V6_n", "V7_p", "V7_n", "V8_p", "V8_n", "V9_p", "V9_n", "V10_p", "V10_n", "GPSString", "HDOP",
                "EXTRANEOUS"]
    outfilename = "{}\\all.txt".format(res_folder)
    # Grab starting and ending points of each survey to reorder
    # NOTE: THIS ASSUMES CONTINUITY WITHIN SURVEY - NO TURNING BOAT AROUND WITHIN SURVEY LINE
    subset = pd.DataFrame(columns=["StartLat", "EndLat", "StartLong", "EndLong", "Filename"])
    excludeSurveys = pd.DataFrame(columns=["Filename", "Number_of_Data_Points"])
    for filename in glob.glob('{}/*.txt'.format(res_folder)):
        if filename == outfilename:
            continue
        # print(filename)
        temp = pd.read_csv(filename, sep=';').reset_index()
        temp.columns = colNames


        # Check if survey is bad (500m or 100 point threshold)
        if len(temp) < 100:
            excludeSurveys = excludeSurveys.append(pd.DataFrame([[filename, str(len(temp))]],
                                                                columns=["Filename",
                                                                         "Number_of_Data_Points"])).reset_index(drop=True)
            logging.info("Resistivity file excluded, length: " + str(len(temp)))
            logging.info(filename + "\n")
            continue

        # Convert the starting and ending coordinates of the survey from degrees decimal minutes to decimal degrees
        try:
            startLat = float(str(temp.loc[0, "Latitude"])[0:2]) + float(str(temp.loc[0, "Latitude"])[2:]) / 60
            endLat = float(str(temp.loc[len(temp)-1, "Latitude"])[0:2]) + float(str(temp.loc[len(temp)-1, "Latitude"])[2:]) / 60
            startLong = float(str(temp.loc[0, "Longitude"])[0:3]) - float(str(temp.loc[0, "Longitude"])[3:]) / 60
            endLong = float(str(temp.loc[len(temp)-1, "Longitude"])[0:3]) - float(str(temp.loc[len(temp)-1, "Longitude"])[3:]) / 60
        except ValueError:
            logging.critical("Could not convert latitude or longitude in " + filename + "\n")
            tkMessageBox.showerror("FORMATTING ERROR",
                                   "Value Error: could not convert latitude or longitude in " + filename)
            exit()
        # Check to see if Longitude is formatted like we want it to
        if startLong > 0 or endLong > 0:
            logging.error("Incorrect longitude format in " + filename + "\n")
            tkMessageBox.showerror("FORMATTING ERROR",
                                   "Error: please format longitude with negative sign for file " + filename)
            exit()
        subset = subset.append(pd.DataFrame([[startLat, endLat, startLong, endLong, filename]],
                                            columns=["StartLat", "EndLat", "StartLong", "EndLong", "Filename"])).reset_index(drop=True)

    # Track files that were removed due to their length
    logging.info("Writing excluded surveys to file\n")
    excludeSurveys.to_csv(path + "\\EXCLUDED_SURVEYS_RES.txt", index=False)

    # Reorganize files based upon their location to one another
    reorderedSubset = pd.DataFrame(columns=["StartLat", "EndLat", "StartLong", "EndLong", "Filename", "Distance", "Reverse"])
    # Pick starting survey as one where start is farthest away from finish
    subset["Distance"] = 0.00
    for i, f in enumerate(subset.Filename):
        startLat = subset.loc[i, "StartLat"]
        startLong = subset.loc[i, "StartLong"]
        dist = 0
        # Find the greatest distance between all lines
        for i2, f2 in enumerate(subset.Filename):
            endLat = subset.loc[i2, "EndLat"]
            endLong = subset.loc[i2, "EndLong"]
            dist = max(dist, haversine(startLong, startLat, endLong, endLat))
        subset.at[i, "Distance"] = dist
    # Start survey is one with greatest starting distance from any survey
    reorderedSubset = reorderedSubset.append(subset.loc[subset["Distance"].idxmax(), :]).reset_index(drop=True)
    reorderedSubset.loc[len(reorderedSubset) - 1, "Reverse"] = False  # First line shouldn't need reversal
    subset.drop([subset["Distance"].idxmax()], inplace=True)
    subset.reset_index(drop=True, inplace=True)

    # Reorder remaining surveys based on distance from the end of previous survey
    while len(subset) > 0:
        endLat = reorderedSubset.loc[len(reorderedSubset)-1, "EndLat"]
        endLong = reorderedSubset.loc[len(reorderedSubset)-1, "EndLong"]
        subset["Distance"] = 999999999.00
        subset["ReverseDistance"] = 999999999.00
        # The next line "starts" closest to the "end" of the previous line
        for i2, f2 in enumerate(subset.Filename):
            startLat = subset.loc[i2, "StartLat"]
            startLong = subset.loc[i2, "StartLong"]
            subset.at[i2, "Distance"] = haversine(startLong, startLat, endLong, endLat)
            subset.at[i2, "ReverseDistance"] = haversine(subset.loc[i2, "EndLong"], subset.loc[i2, "EndLat"], endLong, endLat)

        # Check to see if next survey section is reversed
        if min(subset["Distance"]) <= min(subset["ReverseDistance"]):
            reorderedSubset = reorderedSubset.append(subset.loc[subset["Distance"].idxmin(), :]).reset_index(drop=True)
            reorderedSubset.loc[len(reorderedSubset) - 1, "Reverse"] = False
        else:
            reorderedSubset = reorderedSubset.append(subset.loc[subset["ReverseDistance"].idxmin(), :]).reset_index(drop=True)
            reorderedSubset.loc[len(reorderedSubset)-1, "Reverse"] = True
        subset.drop([subset["Distance"].idxmin()], inplace=True)
        subset.reset_index(drop=True, inplace=True)
    reorderedSubset.loc[0, "Reverse"] = False  # First line shouldn't need reversal (need to restate)

    # Create new filenames for surveys based on their order and export directory to csv file
    reorderedSubset.drop(["StartLat", "EndLat", "StartLong", "EndLong", "Distance", "ReverseDistance"], axis=1, inplace=True)
    reorderedSubset["NewFilename"] = reorderedSubset.index + 1
    reorderedSubset["NewFilename"] = directory + "\\" + userRiverName + "_" + reorderedSubset["NewFilename"].apply(lambda k: str(k).zfill(3)) + ".txt"
    logging.info("Writing renamed resistivity directory to file\n")
    reorderedSubset.to_csv(path + "\\RENAMED_RESISTIVITY_FILE_DIRECTORY.txt", index=False)

    # Rename the actual files in the renamed directory
    logging.info("Copying renamed files to new directory\n")
    for i, fOld in enumerate(reorderedSubset["Filename"]):
        fNew = path + "\\" + reorderedSubset.loc[i, "NewFilename"].replace('/', '\\').split('\\')[-1]
        copyfile(fOld, fNew)

    # %% -----------------------------------------------------------------------------------------------------------------
    # Preprocessing Resistivity Data
    print('Aggregating raw data files')
    logging.info("Aggregating raw data files\n")

    # Copy resistivity data into a single file
    colNames = ["Distance", "Depth", "Rho 1", "Rho 2", "Rho 3", "Rho 4", "Rho 5", "Rho 6", "Rho 7", "Rho 8", "Rho 9",
                "Rho 10", "C1", "C2", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "Latitude",
                "Longitude", "In_p", "In_n", "V1_p", "V1_n", "V2_p", "V2_n", "V3_p", "V3_n", "V4_p", "V4_n", "V5_p", "V5_n",
                "V6_p", "V6_n", "V7_p", "V7_n", "V8_p", "V8_n", "V9_p", "V9_n", "V10_p", "V10_n", "GPSString", "UTC",
                "Latitude2", "D1", "Longitude2", "D2", "Fix Quality", "Satellites", "HDOP", "Altitude", "D3",
                "Height of Geoid", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11"]
    importfile = pd.DataFrame(columns=colNames)
    for i, filename in enumerate(reorderedSubset["Filename"]):
        # If file flagged for reversal, reverse
        if reorderedSubset.loc[i, "Reverse"]:
            temp = pd.read_csv(filename, sep=';|,', engine='python').reset_index()
            temp.columns = colNames
            temp = temp.iloc[::-1]  # Reversal line
            temp["Filename"] = reorderedSubset.loc[i, "NewFilename"]
            bname = os.path.splitext(filename)[0]

            header_size=26
            num_data_bytes=14

            with open('{}.bin'.format(bname),'rb') as fin:
                header = fin.read(header_size)
                data_str = fin.read(num_data_bytes)
                data = data_str.split()[0]

            temp['Date']=data

            temp['Date'] = pd.to_datetime(temp['Date'])
            temp['Date'] = temp['Date'].dt.dayofyear

        # Otherwise, write file to master file normally
        else:
            temp = pd.read_csv(filename, sep=';|,', engine='python').reset_index()
            temp.columns = colNames
            temp["Filename"] = reorderedSubset.loc[i, "NewFilename"]

            bname = os.path.splitext(filename)[0]

            header_size=26
            num_data_bytes=14

            with open('{}.bin'.format(bname),'rb') as fin:
                header = fin.read(header_size)
                data_str = fin.read(num_data_bytes)
                data = data_str.split()[0]

            temp['Date']=data

            temp['Date'] = pd.to_datetime(temp['Date'])
            temp['Date'] = temp['Date'].dt.dayofyear

        # Attempt to remove overlapping lines
        # Creates a box with the corners being the end of the previous line and start of next line
        # Removes any of the next line within that box
        try:
            # Grab the end coordinates of the previous line
            pL_lat = importfile["Latitude"].iloc[-1]
            pL_lon = importfile["Longitude"].iloc[-1]
            # ... And the start coordinates of the current line
            cL_lat = temp.loc[0, "Latitude"]
            cL_lon = temp.loc[0, "Longitude"]
            prevLen = len(temp)
            temp = temp[(temp["Latitude"] > max(pL_lat, cL_lat)) | (temp["Latitude"] < min(pL_lat, cL_lat)) |
                        (temp["Longitude"] > max(pL_lon, cL_lon)) | (temp["Longitude"] < min(pL_lon, cL_lon))]
            if prevLen != len(temp):
                logging.info(str(prevLen-len(temp)) + " overlapping point(s) removed from file " + filename)
        except IndexError:
            # Catches the case where we are looking at the first line - simply skip the removal process
            pass
        importfile = importfile.append(temp)
    importfile.reset_index(drop=True, inplace=True)

    # Write combined file
    importfile.to_csv(outfilename, index=False)

    print('Processing resistivity data')
    logging.info("Processing resistivity data\n")
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
                 'Final_Altitude')
    # Add additional data columns listed above and remove unwanted data columns
    for x in data_cols:
        importfile[x]=np.nan
    importfile.drop(['D1', 'D2', 'D3', "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11"],
                    inplace=True, axis=1)

    importfile.set_index([range(0,len(importfile.Distance))], inplace=True)
    importfile["Latitude"] = importfile["Latitude"].astype(str)
    importfile["Longitude"] = importfile["Longitude"].astype(str)
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
    # Remove erroneous GPS measurements
    importfile = importfile[importfile["Lat"] != 0]
    importfile = importfile[importfile["Lon"] != 0]
    importfile.reset_index(inplace=True, drop=True)

    # Reformat numbers in the file to float
    for col in importfile.columns[1:]:
        try:
            importfile[col] = importfile[col].astype(float)
        except:
            pass

    # %% -----------------------------------------------------------------------------------------------------------------
    # Adding File column to process data in Oasis in chunks
    importfile['File']="1"
    j=1
    for x in range(1,len(importfile.Filename)):
        if importfile.ix[x,"Filename"]==importfile.ix[x-1,"Filename"]:
            importfile.ix[x,"File"]=importfile.ix[x-1,"File"]

        elif importfile.ix[x,"Filename"]!=importfile.ix[x-1,"Filename"]:
            j+=1
            importfile.ix[x,"File"]=j
        else:
            break

    # %% -----------------------------------------------------------------------------------------------------------------
    # Import INI file
    try:
        ini = pd.read_csv(ini_file, index_col=None, sep='=')
        depthoffset = float(ini.ix['DepthOffset', '[SwitchPro]'])
        if depthoffset > 0:
            tkMessageBox.showwarning("WARNING", "Positive value for depth offset from INI file")
            logging.warning("Positive value for depth offset from ini file\n")
    except:
        tkMessageBox.showerror("FILE ERROR", "No INI file selected or incorrect file format...")
        logging.error("No INI file selected or incorrect file format\n")
        exit()

    # %% -----------------------------------------------------------------------------------------------------------------
    # Applying the bandpass filter and rolling average
    logging.info("Applying bandpass filter\n")
    for x in range(1,11):
        band_pass(importfile,'Rho {}'.format(x),0,250)
        rolling_avg(importfile, 'Rho {}'.format(x), 'Rho {}_bandpass'.format(x), 20)

    #%%
    # Applying the depth filter
    logging.info("Applying depth filter\n")
    depth_filt(importfile, 'Depth', depthoffset, 0.01)
    rolling_avg(importfile, 'Depth', 'Depth_filt'.format(x), 20)

#    #%%
#    # Removing large jumps in altitude (commonly caused by bridges)
#    logging.info("Filtering altitude via percentage change\n")
#    importfile['Alt_pct']=importfile['Altitude'].pct_change()
#    importfile['Altitude_bandpass']=importfile['Altitude']
#    importfile['Altitude_bandpass'][importfile['Alt_pct']<np.nanpercentile(importfile['Alt_pct'],5)] = np.nan
#    importfile['Altitude_bandpass'][importfile['Alt_pct']>np.nanpercentile(importfile['Alt_pct'],95)] = np.nan
#    importfile['Altitude_bandpass']=importfile['Altitude_bandpass'].interpolate()


    #%%
    # Filtering Altitude via rolling median filter
    logging.info("Filtering altitude via rolling median filter\n")
    rolling_median(importfile, 'Altitude', 'Altitude', 20)

    #%%
    #Rounding Altitude to the decimeter
    importfile['Altitude_rollmed']=importfile['Altitude_rollmed'].round(1)

    #%%
    # Converting WGS 84 coordinates to UTM 15N coordiantes
    logging.info("Converting WGS84 coordinates\n")
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
    # Calculating the distance from UTM coordinates
    print("Calculating distance from UTM coordinates")
    logging.info("Calculating distance from UTM coordinates\n")
    importfile["Cor_Dist"] = np.sqrt(np.square(importfile['X_UTM'] - importfile['X_UTM'].shift()) +
                                     np.square(importfile['Y_UTM'] - importfile['Y_UTM'].shift()))
#    importfile["Cor_Dist"][0]=0.00
    importfile["Cum_dist"] = importfile["Cor_Dist"].cumsum()
    importfile["Cum_dist"][0]=0.00


    # %% -----------------------------------------------------------------------------------------------------------------
    #Replacing all NaNs with "*" and dropping geometry column
    importfile.fillna('*', inplace=True)
    importfile1 = importfile.drop('geometry', axis=1)

    #%%
    logging.info("Saving processed resistivity file\n")
    saveRes = asksaveasfilename(initialfile='{}_Res.csv'.format(userRiverName),defaultextension='.csv',title="Designate resitivity csv name and location", filetypes=[('csv file', '*.csv')])
    try:
        importfile1.to_csv(saveRes, index=False)
    except IOError:
        logging.critical("Error: could not save resistivity data to file.  Ensure file is not open.")
        tkMessageBox.showerror("FILE ERROR", "Could not save resistivity data to file.  Ensure filename is not open.")
        exit()
    print('Resistivity data exported')


    # %% -----------------------------------------------------------------------------------------------------------------
    # Preprocessing QW Data

    # Grab starting and ending points of each survey to reorder
    # NOTE: THIS ASSUMES CONTINUITY WITHIN SURVEY - NO TURNING BOAT AROUND WITHIN SURVEY LINE
    wqsubset = pd.DataFrame(columns=["StartLat", "EndLat", "StartLong", "EndLong", "Filename"])
    wqexcludeSurveys = pd.DataFrame(columns=["Filename", "Number_of_Data_Points"])
    for filename in glob.glob('{}/*.csv'.format(wq_folder)):
        temp = pd.read_csv(filename, sep=',', skiprows=12, index_col=False, engine='python', encoding='utf-16',
                           names=["Date", "Time", "°C", "mmHg", "DO %", "SPC-uS/cm", "C-uS/cm", "ohm-cm", "pH",
                                  "NH4-N mg/L", "NO3-N mg/L", "Cl mg/L", "FNU", "TSS mg/L", "DEP m", "ALT m", "Lat", "Lon"])

        # Check if survey is bad (only one entry)
        if len(temp) < 2:
            wqexcludeSurveys = wqexcludeSurveys.append(pd.DataFrame([[filename, str(len(temp))]],
                                                                columns=["Filename",
                                                                         "Number_of_Data_Points"])).reset_index(drop=True)
            logging.info("Water quality file excluded, length: " + str(len(temp)))
            logging.info(filename + "\n")
            continue

        # Convert the starting and ending coordinates of the survey from degrees decimal minutes to decimal degrees
        try:
            startLat = temp.loc[0, "Lat"]
            endLat = temp.loc[len(temp)-1, "Lat"]
            startLong = temp.loc[0, "Lon"]
            endLong = temp.loc[len(temp)-1, "Lon"]
        except ValueError:
            logging.critical("Could not convert latitude or longitude in " + filename + "\n")
            tkMessageBox.showerror("FORMATTING ERROR",
                                   "Value Error: could not convert latitude or longitude in " + filename)
            exit()
        # Check to see if Longitude is formatted like we want it to
        if startLong > 0 or endLong > 0:
            logging.error("Incorrect longitude format in " + filename + "\n")
            tkMessageBox.showerror("FORMATTING ERROR",
                                   "Error: please format longitude with negative sign for file " + filename)
            exit()
        wqsubset = wqsubset.append(pd.DataFrame([[startLat, endLat, startLong, endLong, filename]],
                                                columns=["StartLat", "EndLat", "StartLong", "EndLong", "Filename"])).reset_index(drop=True)

    # Track files that were removed due to their length
    logging.info("Writing excluded surveys to file\n")
    wqexcludeSurveys.to_csv(path + r"\\EXCLUDED_SURVEYS_WQ.txt", index=False)

    # Reorganize files based upon their location to one another
    wqreorderedSubset = pd.DataFrame(columns=["StartLat", "EndLat", "StartLong", "EndLong", "Filename", "Distance", "Reverse"])
    # Pick starting survey as one where start is closest to the resistivity start
    wqsubset["Distance"] = 0.00
    startLat = importfile1.loc[0, "Lat"]
    startLong = importfile1.loc[0, "Lon"]
    for i, f in enumerate(wqsubset.Filename):
        endLat = wqsubset.loc[i, "EndLat"]
        endLong = wqsubset.loc[i, "EndLong"]
        dist = haversine(startLong, startLat, endLong, endLat)
        wqsubset.at[i, "Distance"] = dist

    # Start survey is one with shortest starting distance from the first resistivity survey
    wqreorderedSubset = wqreorderedSubset.append(wqsubset.loc[wqsubset["Distance"].idxmin(), :]).reset_index(drop=True)
    wqreorderedSubset["Reverse"] = False  # First line shouldn't need reversal
    wqsubset.drop([wqsubset["Distance"].idxmin()], inplace=True)
    wqsubset.reset_index(drop=True, inplace=True)

    # Reorder remaining surveys based on distance from the end of previous survey
    while len(wqsubset) > 0:
        endLat = wqreorderedSubset.loc[len(wqreorderedSubset)-1, "EndLat"]
        endLong = wqreorderedSubset.loc[len(wqreorderedSubset)-1, "EndLong"]
        wqsubset["Distance"] = 999999999.00
        wqsubset["ReverseDistance"] = 999999999.00
        # The next line "starts" closest to the "end" of the previous line
        for i2, f2 in enumerate(wqsubset.Filename):
            startLat = wqsubset.loc[i2, "StartLat"]
            startLong = wqsubset.loc[i2, "StartLong"]
            wqsubset.at[i2, "Distance"] = haversine(startLong, startLat, endLong, endLat)
            wqsubset.at[i2, "ReverseDistance"] = haversine(wqsubset.loc[i2, "EndLong"], wqsubset.loc[i2, "EndLat"], endLong, endLat)
        # Check to see if next survey section is reversed
        if min(wqsubset["Distance"]) <= min(wqsubset["ReverseDistance"]):
            wqreorderedSubset = wqreorderedSubset.append(wqsubset.loc[wqsubset["Distance"].idxmin(), :]).reset_index(drop=True)
            wqreorderedSubset.loc[len(wqreorderedSubset) - 1, "Reverse"] = False
        else:
            wqreorderedSubset = wqreorderedSubset.append(wqsubset.loc[wqsubset["ReverseDistance"].idxmin(), :]).reset_index(drop=True)
            wqreorderedSubset.loc[len(wqreorderedSubset)-1, "Reverse"] = True
        wqsubset.drop([wqsubset["Distance"].idxmin()], inplace=True)
        wqsubset.reset_index(drop=True, inplace=True)
    wqreorderedSubset.loc[0, "Reverse"] = False  # First line shouldn't need reversal (need to restate)

    # Create new filenames for surveys based on their order and export directory to csv file
    wqreorderedSubset.drop(["StartLat", "EndLat", "StartLong", "EndLong", "Distance", "ReverseDistance"], axis=1, inplace=True)
    wqreorderedSubset["NewFilename"] = wqreorderedSubset.index + 1
    wqreorderedSubset["NewFilename"] = directory + "\\" + userRiverName + "_" + wqreorderedSubset["NewFilename"].apply(lambda k: str(k).zfill(3)) + "_WQ.csv"
    logging.info("Writing renamed water quality directory to file\n")
    wqreorderedSubset.to_csv(path + r"\\RENAMED_WQ_FILE_DIRECTORY.txt", index=False)

    # Rename the actual files in the renamed directory
    logging.info("Copying renamed files to new directory\n")
    for i, fOld in enumerate(wqreorderedSubset["Filename"]):
        fNew = path + "\\" + wqreorderedSubset.loc[i, "NewFilename"].replace('/', '\\').split('\\')[-1]
        copyfile(fOld, fNew)

    # %% -----------------------------------------------------------------------------------------------------------------
    # Import the water quality data based on the reordered index
    print('Importing water quality data')
    wqCols = ["Date", "Time", "°C", "mmHg", "DO %", "SPC-uS/cm", "C-uS/cm", "ohm-cm", "pH", "NH4-N mg/L", "NO3-N mg/L",
              "Cl mg/L", "FNU", "TSS mg/L", "DEP m", "ALT m", "Lat", "Lon"]
    qwdata = pd.DataFrame(columns=wqCols)
    for i, filename in enumerate(wqreorderedSubset["Filename"]):
        # If file flagged for reversal, reverse
        if wqreorderedSubset.loc[i, "Reverse"]:
            temp = pd.read_csv(filename, sep=',', skiprows=12, index_col=False, engine='python', encoding='utf-16',
                               names=wqCols)
            temp = temp.iloc[::-1]  # Reversal line
            temp["Filename"] = wqreorderedSubset.loc[i, "NewFilename"]
        # Otherwise, write file to master file normally
        else:
            temp = pd.read_csv(filename, sep=',', skiprows=12, index_col=False, engine='python', encoding='utf-16',
                               names=wqCols)
            temp["Filename"] = wqreorderedSubset.loc[i, "NewFilename"]
        # Attempt to remove overlapping lines
        # Creates a box with the corners being the end of the previous line and start of next line
        # Removes any of the next line within that box
        try:
            # Grab the end coordinates of the previous line
            pL_lat = qwdata["Lat"].iloc[-1]
            pL_lon = qwdata["Lon"].iloc[-1]
            # ... And the start coordinates of the current line
            cL_lat = temp.loc[0, "Lat"]
            cL_lon = temp.loc[0, "Lon"]
            prevLen = len(temp)
            temp = temp[(temp["Lat"] > max(pL_lat, cL_lat)) | (temp["Lat"] < min(pL_lat, cL_lat)) |
                        (temp["Lon"] > max(pL_lon, cL_lon)) | (temp["Lon"] < min(pL_lon, cL_lon))]
            if prevLen != len(temp):
                logging.info(str(prevLen-len(temp)) + " overlapping point(s) removed from file " + filename)
        except IndexError:
            # Catches the case where we are looking at the first line - simply skip the removal process
            pass
        qwdata = qwdata.append(temp)
    qwdata.reset_index(drop=True, inplace=True)
    print('Water-quality data imported')

    # %% -----------------------------------------------------------------------------------------------------------------
    print('Processing water-quality data')
    logging.info("Processing water-quality data\n")
    # Remove erroneous GPS measurements
    qwdata = qwdata[qwdata["Lat"] != 0.00]
    qwdata = qwdata[qwdata["Lon"] != 0.00]
    qwdata.reset_index(inplace=True, drop=True)

    qwdata.dropna(axis=1, how='all', inplace=True)
    qwdata.rename(columns={'°C':'Temp_C', 'SPC-uS/cm': 'SPC_mscm','C-uS/cm':'Cond_mscm','ohm-cm':'Res_ocm','ALT m':'Alt_m'}, inplace=True)
    try:
        qwdata = qwdata.loc[qwdata.Res_ocm!="    +++++",:]
    except:
        pass
    qwdata['Res_ocm'] = qwdata['Res_ocm'].astype('float')
    qwdata['Ohm_m']=qwdata['Res_ocm']/100

    # %% -----------------------------------------------------------------------------------------------------------------
    # Applying a rolling average on resistivity
    rolling_avg(qwdata, 'Ohm_m', 'Ohm_m', 20)

    #%%
    for col in qwdata.columns[1:]:
        try:
            qwdata[col] = qwdata[col].astype(float)
        except:
            pass

    # %% -----------------------------------------------------------------------------------------------------------------
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

    # %% -----------------------------------------------------------------------------------------------------------------
    # Adding File column to process data in Oasis in chunks
    qwdata['File']="1"
    j=1
    for x in range(1,len(qwdata.Filename)):
        if qwdata.ix[x,"Filename"]==qwdata.ix[x-1,"Filename"]:
            qwdata.ix[x,"File"]=qwdata.ix[x-1,"File"]

        elif qwdata.ix[x,"Filename"]!=qwdata.ix[x-1,"Filename"]:
            j+=1
            qwdata.ix[x,"File"]=j
        else:
            break

    #%%
    #Exporting resistivity data as a shapefile
    logging.info("Saving processed water-quality shapefile\n")
    saveWQshp = asksaveasfilename(initialfile='{}_WQ.shp'.format(userRiverName),defaultextension='.shp',title="Designate water-quality shapefile name and location", filetypes=[('shp file', '*.shp')], initialdir=directory)
    try:
        qwdata.to_file(saveWQshp,driver='ESRI Shapefile')
    except IOError:
        logging.critical("Error: could not save processed water-quality data to shapefile.  Ensure file is not open.")
        tkMessageBox.showerror("FILE ERROR", "Could not save processed water-quality data to shapefile.  Ensure filename is not open.")
        exit()
    print('Processed water-quality shapefile exported')

    #%%
    # Creating buffers and spatially joining QW with resitivity data
    Ohm_buffer = qwdata.buffer(5)
    qwdata1 = qwdata[['X_UTM','Y_UTM','Ohm_m_rollavg','Temp_C','Date','Time','geometry']]
    qwdata1['geometry'] = Ohm_buffer
    qwdata1 = qwdata1.set_geometry('geometry')
    resOhm = gp.sjoin(importfile,qwdata1,how='left', op='intersects')
    resOhm[['Ohm_m_rollavg','Temp_C']] = resOhm[['Ohm_m_rollavg','Temp_C']].interpolate()
    resOhm[['Ohm_m_rollavg','Temp_C']] = resOhm[['Ohm_m_rollavg','Temp_C']].fillna(method='bfill')
    resOhm['Temp_C'] = resOhm['Temp_C'].round(1)

    #%%
    # Populating the final fields
    resOhm['Ohm_m'] = resOhm['Ohm_m_rollavg']
    resOhm['Final_Altitude'] = resOhm['Altitude_rollmed']
    resOhm['Cor_Depth'] = resOhm['Depth_filt']
    resOhm['Final_Rho_1'] = resOhm['Rho 1_rollavg']
    resOhm['Final_Rho_2'] = resOhm['Rho 2_rollavg']
    resOhm['Final_Rho_3'] = resOhm['Rho 3_rollavg']
    resOhm['Final_Rho_4'] = resOhm['Rho 4_rollavg']
    resOhm['Final_Rho_5'] = resOhm['Rho 5_rollavg']
    resOhm['Final_Rho_6'] = resOhm['Rho 6_rollavg']
    resOhm['Final_Rho_7'] = resOhm['Rho 7_rollavg']
    resOhm['Final_Rho_8'] = resOhm['Rho 8_rollavg']
    resOhm['Final_Rho_9'] = resOhm['Rho 9_rollavg']
    resOhm['Final_Rho_10'] = resOhm['Rho 10_rollavg']

    #%%
    # Converting *s to NaNs for import into GIS as a float
    resOhm.replace('*',np.nan, inplace=True)

    #%%
    logging.info("Saving preliminary merged QW/resistivity shapefile\n")
    savepreres = asksaveasfilename(initialfile='{}_Merged_QWRes.shp'.format(userRiverName),defaultextension='.shp',title="Designate preliminary merged QW/resitivity shapefile name and location", filetypes=[('shp file', '*.shp')], initialdir=directory)
    try:
        resOhm.to_file(savepreres,driver='ESRI Shapefile')
    except IOError:
        logging.critical("Error: could not save preliminary merged QW/resistivity data to shapefile.  Ensure file is not open.")
        tkMessageBox.showerror("FILE ERROR", "Could not save preliminary merged QW/resistivity data to shapefile.  Ensure filename is not open.")
        exit()
    print('Preliminary merged QW/resistivity shapefile exported')

    #%%
    resOhm_df = resOhm.drop('geometry',axis=1)

    #%%
    logging.info("Export preliminary merged QW/resisitivty data\n")
    resOhm_csv = asksaveasfilename(initialfile='{}_Merged_WQRes.csv'.format(userRiverName),defaultextension='.csv',title="Designate preliminary merged QW/resisitivty csv name and location", filetypes=[('csv file', '*.csv')], initialdir=directory)
    try:
        resOhm_df.to_csv(resOhm_csv, index=False)
    except IOError:
        logging.critical("Error: could not save preliminary merged QW/resisitivty data to file.  Ensure file is not open.")
        tkMessageBox.showerror("FILE ERROR", "Could not save preliminary merged QW/resisitivty data to file.  Ensure filename is not open.")
        exit()
    print('Preliminary merged QW/resistivity csv exported!')

    # %% -----------------------------------------------------------------------------------------------------------------
    # Removing geometry to dataframe
    qwdata.drop('geometry',axis=1,inplace=True)

    #%%
    logging.info("Export water quality data\n")
    saveQW = asksaveasfilename(initialfile='{}_WQ.csv'.format(userRiverName),defaultextension='.csv',title="Designate water-quality csv name and location", filetypes=[('csv file', '*.csv')], initialdir=directory)
    try:
        qwdata.to_csv(saveQW, index=False)
    except IOError:
        logging.critical("Error: could not save water-quality data to file.  Ensure file is not open.")
        tkMessageBox.showerror("FILE ERROR", "Could not save water-quality data to file.  Ensure filename is not open.")
        exit()
    print('Water-quality csv exported!')

    # %% -----------------------------------------------------------------------------------------------------------------
    # Write summary file
    try:
        summaryFile = open(directory + "\\OASIS_PREPROCESSING_SUMMARY.txt", "w+")
    except IOError:
        logging.critical("Error: could not write summary file")
        exit()
    summaryFile.write("Processed on {:%Y-%m-%d %H:%M:%S}\n\n".format(datetime.datetime.now()))
    totalDistance = importfile1.loc[len(importfile1)-1, "Cum_dist"]/1000
    summaryFile.write("Total distance processed: %.2f" % totalDistance + " kilometers\n")
    summaryFile.write("Number of resistivity files read: " + str(len(reorderedSubset)) + "\n")
    for f in reorderedSubset.NewFilename:
        summaryFile.write(f)
        summaryFile.write('\n')
    summaryFile.write("\n\nNumber of water quality files read: " + str(len(wqreorderedSubset)) + "\n")
    for f in wqreorderedSubset.NewFilename:
        summaryFile.write(f)
        summaryFile.write('\n')
    summaryFile.write("\n\n")
    summaryFile.close()
    #%%
    # Data Release
    if eg.ynbox(title='Data Release Utility',msg='Do you want to export raw and processed data in a data release format?'):
        fieldNames=['Iris Serial Number','Cable Serial Number','EchoSounder GPS Serial Number','QW Probe Serial Number']
        msg='Fill in the values'
        title='Data Release Form'
        fieldValues=['18705-4177458684-507','0117352673-87','1532702SCSC','16F102579']
        eg.multenterbox(msg,title,fieldNames, fieldValues)
        if fieldValues is None:
            sys.exit(0)

        while 1:
            errmsg = ""
            for i, name in enumerate(fieldNames):
                if fieldValues[i].strip() == "":
                    errmsg += "{} is a required field.\n\n".format(name)
            if errmsg == "":
                break # no problems found

            fieldValues = eg.multenterbox(errmsg, title, fieldNames, fieldValues)
            if fieldValues is None:
                break

        #%%
        dr_raw=importfile[['File','Date','UTC','Depth','Lat','Lon','Altitude','Cum_dist','In_n','In_p','V1_n','V1_p','V2_n','V2_p','V3_n','V3_p','V4_n','V4_p','V5_n','V5_p','V6_n','V6_p','V7_n','V7_p','V8_n','V8_p','V9_n','V9_p','V10_n','V10_p','Rho 1','Rho 2','Rho 3','Rho 4','Rho 5','Rho 6','Rho 7','Rho 8','Rho 9','Rho 10','C1','C2','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11']]

        dr_post=importfile[['File','Date','UTC','Depth_rollavg','Ohm_m','Lat','Lon','Altitude_rollmed','Cum_dist','Rho 1_rollavg','Rho 2_rollavg','Rho 3_rollavg','Rho 4_rollavg','Rho 5_rollavg','Rho 6_rollavg','Rho 7_rollavg','Rho 8_rollavg','Rho 9_rollavg','Rho 10_rollavg']]
        #%%
        dr_raw['Iris_SN']=fieldValues[0]
        dr_raw['Cable_SN']=fieldValues[1]
        dr_raw['Echo_GPS_SN']=fieldValues[2]
        dr_raw['QW_SN']=fieldValues[3]

        dr_post['Iris_SN']=fieldValues[0]
        dr_post['Cable_SN']=fieldValues[1]
        dr_post['Echo_GPS_SN']=fieldValues[2]
        dr_post['QW_SN']=fieldValues[3]

        dr_raw.rename(columns={'File':'Profile','UTC':'Time','Lat':'Latitude','Lon':'Longitude','Cum_dist':'UTM_distance','Rho 1':'Rho_1','Rho 2':'Rho_2','Rho 3':'Rho_3','Rho 4':'Rho_4','Rho 5':'Rho_5','Rho 6':'Rho_6','Rho 7':'Rho_7','Rho 8':'Rho_8','Rho 9':'Rho_9','Rho 10':'Rho_10','Altitude':'Elevation'}, inplace=True)

        dr_post.rename(columns={'File':'Profile','UTC':'Time','Lat':'Latitude','Lon':'Longitude','Cum_dist':'UTM_distance','Rho 1_rollavg':'Rho1','Rho 2_rollavg':'Rho2','Rho 3_rollavg':'Rho3','Rho 4_rollavg':'Rho4','Rho 5_rollavg':'Rho5','Rho 6_rollavg':'Rho6','Rho 7_rollavg':'Rho7','Rho 8_rollavg':'Rho8','Rho 9_rollavg':'Rho9','Rho 10_rollavg':'Rho10','Altitude':'Elevation','Ohm_m':'Water_Res'}, inplace=True)

        raw = eg.filesavebox(title="Save raw data release file as...",default='{}_Raw_DataRelease.csv'.format(userRiverName),filetypes=['*.csv'])
        post = eg.filesavebox(title="Save prcoessed data release file as...",default='{}_Processed_DataRelease.csv'.format(userRiverName),filetypes=['*.csv'])
        dr_post.to_csv(post, index=False)
        dr_raw.to_csv(raw, index=False)
    else:
       if choice=="Oasis Preprocessor":
           sys.exit(0)
       else:
           pass

#%%
#Bring up GUI and execute functions
def catchEmAll(*exc_info):
    errormessage = "".join(traceback.format_exception(*exc_info))
    logging.critical("Uncaught error encountered: \n%s", errormessage)

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
# Record logging events to log file
logging.basicConfig(filename="\\OASIS_PREPROCESSING_LOGFILE.txt", format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p', filemode='w', level=logging.INFO)

sys.excepthook = catchEmAll

ret_val = eg.msgbox("USGS Workbench Preprocessor and Data Release Utility")
if ret_val is None: # User closed eg.msgbox
    sys.exit(0)

title ="USGS Oasis/Workbench Preprocessor and Data Release Utility"
msg = "Choose which utility you would like to use"
choices = ["Oasis Preprocessor","Workbench Preprocessor","Oasis/Workbench Preprocessor"]

while 1:
    choice = eg.buttonbox(msg, title, choices)
    if choice is None:
        sys.exit(0)
    elif choice=="Oasis Preprocessor":
        oasis()
    elif choice=="Workbench Preprocessor":
        workbench()
        workbench_checks()
    elif choice=="Oasis/Workbench Preprocessor":
        oasis()
        workbench()
        workbench_checks()

