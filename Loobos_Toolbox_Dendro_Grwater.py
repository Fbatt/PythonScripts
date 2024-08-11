#--- import packages
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import plotly.express as px
#import cufflinks as cf
from datetime import datetime, timedelta
import re

#--- define functions to match regex
def notmatch(infilename):
    raise Exception('while processing %s, encountered columns with names not matching pattern "jrXXXX".'%(infilename))

#--- define functions to import dendrometer data 
def Read_Dendrometers(datapath):
    df_dendro      = pd.DataFrame() #final container object for all years of all dendrometers
    for iname in ['I5','K7','K8','K9']: #dendrometer O3 is not included because the dates match a different pattern
        print('Loading dendrometer %s'%iname)
        infilename  = '%sLoobos_dendrometers\%s.csv'%(datapath,iname) # Profile data without storage or gapfilling
        print('Loading %s'%infilename)
        df_tmp = pd.read_csv(infilename,sep=",",index_col=False)

        #the expected output of read_csv is a dataframe with in the first column groeiseizoen and 
        #testing the assumption that the first column is labelled 'groeiseizoen'
        if df_tmp.columns[0]!='groeiseizoen':
            raise Exception('first column of %s is not labelled "groeiseizoen".'%(infilename))
        #test to see if all column names beside first match the pattern 'jrXXXX'
        colnames = df_tmp.columns[1:]
        for colname in colnames:
            re.match(r"^jr[0-9]{4}$", colname) or notmatch(infilename)  #this works because the parser evaluates re.match and if it can't (doesn't match pattern), only then does it evaluate what's after the or.
        #trim the 'jr' part off the column names
        years = [colname[2:] for colname in colnames]

        #initialize the dataframe for one year:
        df_dendroyears=pd.DataFrame() #temporary container object for one year of one dendrometer 
        df_year=[] #temporary container
        year_date_time_series= years[0] + "-" + df_tmp['groeiseizoen']
        year_date_time_series=pd.to_datetime(year_date_time_series, format= '%Y-%m-%d %H:%M')
        #df_year = pd.concat([year_date_time_series, df_tmp[colnames[0]]], axis=1)
        #df_year = pd.DataFrame({'dendro_growth': df_tmp[colnames[0]]},index=pd.Index(year_date_time_series, name='datetime')) #dunno why this doesn't work so doing it the difficult way below:
        df_year = pd.DataFrame({'datetime': year_date_time_series ,'dendro_growth_'+iname: df_tmp[colnames[1]]})
        df_year.set_index( df_year['datetime'], inplace=True)
        df_year=df_year.drop(['datetime'], axis=1)
        df_dendroyears=df_year #set up first year of dendroyears so we can join things to it.
        
        #now that we have a start, loop over the other years and vertically concatenate
        for year in years[1:]:
            df_year=[] #temporary container
            print('loading year %s'%(year))
            columntag='jr'+year #produces something like 'jr2008'
            year_date_time_series= year + "-" + df_tmp['groeiseizoen']
            year_date_time_series=pd.to_datetime(year_date_time_series, format= '%Y-%m-%d %H:%M')
            df_year = pd.DataFrame({'datetime': year_date_time_series ,'dendro_growth_'+iname: df_tmp[columntag]}) #we are targetting columntag instead of first instance of colnames like above
            df_year.set_index( df_year['datetime'], inplace=True)
            df_year=df_year.drop(['datetime'], axis=1)
            #df_year is now ready to be joined with df_dendroyears vertically
            df_dendroyears=pd.concat([df_dendroyears,df_year],axis=0) #axis=0 specifies concat on index.
        
        #we join this information to df_dendro and continue iterating over all the dendrometer files
        if df_dendro.empty == True: #if this is the first iteration, copy dendroyears
            df_dendro=df_dendroyears.copy()
            #print('df is empty!')
        elif df_dendro.empty == False: #if it isnt, then add corresponding data to df_dendro
            df_dendro = df_dendro.merge(df_dendroyears, how='outer', left_index = True, right_index = True)
        else:
            raise Exception('something went wrong with appending a single year to df_dendro, consult code to fix')
        #print(df_dendro)
        
    return df_dendro

#--- define functions to import groundwater
#this part is still work in progress


