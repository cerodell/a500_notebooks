#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 13:47:36 2017

This will retreive a sounding from U of Wyoming's Web Server
and save it to a text file. This was coded to work with SHARPY

This program will cut out the header, below ground data, and the extra indices 
found after the sounding retrievals making a plain text file with columns. 
This will fill empty data with -9999

You will find two text files are created. The first is the raw data from Wyoming's
server. The second, with suffix: sharppy.txt is the one to import into SHARPpy.
You really don't need the raw data except maybe the data at the end of the file
might be of interest.

@author: Robinson Wallace
@contact: wallacer@colorado.edu
"""
from bs4 import BeautifulSoup
import requests
import pandas as pd
import os
import errno



user      = 'rodell'

#filein  = '/Users/'+user+'/Google Drive File Stream/Shared drives/Research/CRodell/BL/Data/'
#save    = '/Users/'+user+'/Google Drive File Stream/Shared drives/Research/CRodell/BL/Images/'



def make_sure_path_exists(path):
    
    """
    Name: make_sure_path_exists
    
    Inputs: (String) A directory path 
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def station_number(soundingid):
    
    """
    Name: station_num
    Purpose: To return the station number associated with the 4 letter code
             of the station ID used by U of Wyoming's radiosonde archive
    Inputs:
        soundingid: (String) 4 letter code of the radiosonde station
        
    Notes: This list needs to be expanded to include the stations the user
           wants to review. Bobby was too lazy to add them all, so he leaves
           it up to you to add the stations you're interested in. 
    """
    
    
    if soundingid == "KDNR":
        soundingnum = '72469'
#    if soundingid == "KRAP":
#        soundingnum = '72662'
#    if soundingid == "KLBF":
#        soundingnum = '72562'
#    if soundingid == "ABQ":
#        soundingnum = '72365'
#    if soundingid == "FGZ":
#        soundingnum = '72376'
        
        
    return soundingnum

def write_UofWweb(year,month,day,hour,soundingid,save_directory):
    
    """
    Name: write_UofWweb
    
    Purpose: 1) To retireve operational radiosonde data from the University of 
                Wyoming's data archive
             2) Write the sounding in a format SHARPpy can ingest.
    
    Notes: 1) This will leave the retrieved, unformatted data, in the save
              directory. The correctly formatted file is tagged with "_sharppy"
              in the filename. 
           
           2) This should probably get broken up into two functions so the 
              user can just get the radiosonde data, but I think the purpose
              of this routine is to get the data in a format SHARPpy can use,
              so I kept this "simple" by reducing the number of functions.
              
    Areas of Possible Improvement: 
        
            1) This currently only gets soundings that were launched at 
               00 or 12 UTC. Some soundings are indeed available at other hours.
               Adding this functionality is an area of possible improvement.
              
            2) A call to a date that doesn't exist (e.g. August 32) will
               just return an nearly empty text file. A condition could be
               written in to let the user know they input an invalid date.
    """
    
    
    soundingnum = station_number(soundingid)
    
    year=str(year)    
    month = str(month).zfill(2)
    day=str(day).zfill(2)
    
    #Make sure the user put in an hour that exists in Wyoming's database
    if float(hour) != 0 or float(hour) !=12:
        if float(hour) < 12:
            Zhour = '00'
        else:
            Zhour = '12'

    # Assign the web address to UofWyoming's Sounding based on input
    link="http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&"\
         "YEAR="+year+"&MONTH="+month+"&FROM="+day+Zhour+"&TO="+day+Zhour+"&STNM="+soundingnum

    soup = BeautifulSoup(requests.get(link).text, "lxml")
    lines = soup.get_text().split('\n')
        
    fname = soundingid+"_"+year+month+day+"_"+Zhour.zfill(2)+"UTC.txt"
    
    # Make the save directory if it doesn't exist
    make_sure_path_exists(save_directory)
    
    with open(save_directory+fname, 'w', newline='') as f:
        for line in lines:
            line = line+'\n'
            f.write(line) 

    names=['PRES','HGHT','TEMP','DWPT','RH','W','WDIR','WSPD','THETA','THTE','THETV']
    df0 = pd.read_table(save_directory+fname,skiprows=10,skipfooter=54,engine='python',
                       names=names, header=None, sep='\s+',index_col=False)
    df0.fillna(-9999, inplace=True)
    
#    hedr = ' LEVEL \t HGHT \t TEMP \t DWPT \t WDIR \t WSPD'
#    df0['WDIR'][df0['WDIR']==360] = 0
#    df0['WDIR'][df0['WDIR']>360]  = -9999
#    df0['WSPD'][df0['WSPD']<0]    = -9999
#    
#    # Format the columns correctly
#    df1 = df0[['PRES','HGHT','TEMP','DWPT','WDIR','WSPD']]
#    for key in df1:
#        if key == 'WSPD':
#            df1[key] = df1[key].astype(str)+'0'
#            break
#        if key == 'HGHT':
#            df1[key] = df1[key].astype(str)+'.00,'
#        else:
#            df1[key] = df1[key].astype(str)+'0,'
#    
#    # Make sure heights are sequentially higher. If not, replace with -9999
#    for h in range(len(df1['HGHT'])-1):
#        if float(float(df1["HGHT"][h+1].split(',')[0])) <= float(float(df1["HGHT"][h].split(',')[0])):
#            df1['HGHT'][h] = '-9999,'
#    
#    fname = soundingid+"_"+year+month+day+"_"+Zhour.zfill(2)+"UTC_sharppy.txt"
#    with open(save_directory+fname, 'w') as fo:
#        fo.write('%TITLE%\n')
#        fo.write(' '+soundingid+'   '+year[2:4]+month+day+'/'+hour.zfill(4)+'\n')
#        fo.write('\n')
#        fo.write(hedr+'\n')
#        fo.write('-------------------------------------------------------------------\n')
#        fo.write('%RAW%\n')
#        with pd.option_context('display.max_rows', 200, 'display.max_columns', 200):
#            fo.write(df1.to_string(index=False,header=False,justify='right'))
#        fo.write('\n%END%')
#        

if __name__ == "__main__":

    # Format for each entry is YYYYMMDD_HH_SoundingStationID
    # Sounding station ID is the four letter code for the station used by 
    # Wyoming's sounding archive. For example, Denver is coded as KNDR.
    # Use 00 or 12 for the hour. Any other input hour will be rounded to 
    # either 00 or 12.
    # An example for this variable is: datelist = ['20160831_00_KDNR']
    # Get codes from http://weather.uwyo.edu/upperair/sounding.html
    datelist = ['20160630_00_KDNR']
    
    # Identify the directory to save the soundings
    
    
    save_directory = '/Users/'+user+'/Google Drive File Stream/Shared drives/Research/CRodell/BL/Data/'

    
    for d in range(len(datelist)):
        year = datelist[d][0:4].zfill(4)
        month = datelist[d][4:6].zfill(2)
        day = datelist[d][6:8].zfill(2)
        hour = datelist[d][9:11].zfill(2)
        soundingid = datelist[d][12:]
        write_UofWweb(year,month,day,hour,soundingid,save_directory)