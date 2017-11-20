
# coding: utf-8

# In[31]:

"""
Last revised 09/01/2017

Jordan Wilson

This script takes the resistivity survey information from Oasis Montaj and formats it for import into Workbench.

"""

# In[2]:
from Tkinter import *
import pandas as pd
from Tkinter import Tk
from tkFileDialog import asksaveasfilename
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename
import tkMessageBox


#%%
# Insert tkinter buttons here to fork the script between processing for data release output or workbench output. 
top = Tkinter.Tk()
def datarelease():
    print('in progress')
        
def workbench():
    root = Tk() # we don't want a full GUI, so keep the root window from appearing
    root.lift()
    tkMessageBox.showinfo("Directions", "For this script to work, the final Rho columns must be named 'Final_Rho_1,Final_Rho_2,...,Final_Rho_n', Latitude and Longitude columns must be named 'Lat' and 'Lon', and the corrected distance and depth columns should be named 'Cor_Dist' and 'Cor_Depth', and the QW columns containing resisitivty values from the QW meter should be named 'Ohm_m'")
    infile = askopenfilename(title="Select Oasis output csv file for processing") # show an "Open" dialog box and return the path to the selected file
    Tk().withdraw
        
    try:
        data = pd.read_csv(infile)
        print('Oasis file read')
    except:
        print('Oasis file not read')
        
    try:
        data['File_w']=data['Line'].str.split('L',1)
        data['File']=''
        data['File']=data.apply(lambda row: row['File_w'][1],axis=1)
        data['File']=data['File'].astype('int')
    
        alt_col = [col for col in data.columns if 'Alt' in col]
    
        out = pd.DataFrame()
        
    except KeyError,e:
        tkMessageBox.showerror("Missing column",'The following column is missing in the input file: %s. Check to make sure all required columns are present.' % str(e))
    
    data.rename(columns={'File':'Profile'}, inplace=True)

    headers=('/Water_Res','Cor_Dist','Cor_Depth','Rho_1','Rho_2','Rho_3','Rho_4','Rho_5','Rho_6','Rho_7','Rho_8','Rho_9','Rho_10','C1','C2','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','Lat','Long','Altitude','Profile')
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
                 'Long',
                 'Final_Altitude',
                 'Profile')
    
    zipped = zip(headers, data_cols)
    try:
        for x,y in reversed(zipped):
            out.insert(0,x,data['{}'.format(y)].values)
        print('File formatted')
        outfile = asksaveasfilename(title="Save processed file as...",filetypes=[('csv file', '*.csv')])
        out.to_csv(outfile,index=False)
        print('Output csv file saved')
    except KeyError, e:
       tkMessageBox.showerror("Missing column",'The following column is missing in the input file: %s. Check to make sure all required columns are present.' % str(e))
        
    root.destroy()
B1 = Tkinter.Button(top, text = "Data Release", command = datarelease)
B2 = Tkinter.Button(top, text = "Workbench", command = workbench)
L1 = Tkinter.

B1.pack()
B2.pack()
top.mainloop()
# In[2]:

raw_input('Press ENTER to exit')


# In[ ]:



