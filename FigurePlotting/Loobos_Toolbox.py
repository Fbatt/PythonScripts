#--- import packages
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import plotly.express as px
#import cufflinks as cf
from datetime import datetime, timedelta

#--- define functions to construct time array
def dateparse(x):
    year    = int(x[  : 4])
    doy     = int(x[ 4: 7])
    hours   = int(x[ 7: 9])
    minutes = int(x[ 9:11])
    loobos_date = datetime(year, 1, 1) + timedelta(days=doy-1, hours=hours, minutes=minutes)
    return loobos_date

def dateparse_Gapfilled(x):
    year    = int(x[  : 4])
    month   = int(x[ 5: 7])
    day     = int(x[ 8:10])
    hours   = int(x[11:13])
    minutes = int(x[14:16])
    seconds = int(x[17:19])
    loobos_date = datetime(year, month, day) + timedelta(hours=hours, minutes=minutes)
    return loobos_date    
    
#def dateparse_Comb(year,doy,hour):
#    print(year)
#    print(doy)
#    print(hour)
#    YearDayTime = pd.to_datetime(year
   #year    = int(year)
   #doy     = int(doy)
#    hh      = np.floor(float(hour))
#    mm      = np.mod  (float(hour),1)*60
#    
#    if doy > 366:
#        doy  = doy-365
#        year = year + 1
#    timestr = '%4d-%03d-%02d:%02d'%(year,doy,hh,mm)
#
#    YearDayTime = pd.to_datetime(timestr,format='%Y-%j-%H:%M')
##   YearDayTime = pd.to_datetime(timestr).strftime('%Y-%j-%H:%M')
#   return #YearDayTime    


#--- 1. Load raw Eddy Covariance data (without storage and without gapfilling)
def Read_LoobosEddFinal(years,datapath):
    df_EC           = []
    for iy in years:
        print('Loading %d'%iy)
        infilename  = '%s\Eddy\Final\LoobosEdd%02dfinal.csv'%(datapath,np.mod(iy,100)) # EC data without storage or gapfilling
       #df_EC         .append(pd.read_csv(infilename, index_col='YearDayTime', parse_dates=['YearDayTime'], date_parser=dateparse))
       #df_EC         .append(pd.read_csv(infilename, index_col='YearDayTime', parse_dates=['YearDayTime'],date_format='%Y%m%d%H%M'))
        df_EC         .append(pd.read_csv(infilename, index_col='YearDayTime'))
    df_EC           = pd.concat(df_EC,   axis=0, ignore_index=False)
    
    year            = np.asarray(np.floor(df_EC.index/10000000).astype(int)-1970, dtype='datetime64[Y]')
    doy             = np.asarray(         df_EC['Doy'  ]                        , dtype='timedelta64[D]')
    hour            = np.asarray(np.floor(df_EC['Dtime']/100)                   , dtype='timedelta64[h]')
    mins            = np.asarray(np.mod(  df_EC['Dtime']/100,1)*100             , dtype='timedelta64[m]')
    df_EC.index     = pd.to_datetime(year+doy-1+hour+mins)
    df_EC.index.name = 'YearDayTime'
    
    # Make sure the time index is complete, so all data frames have the same length.
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_EC           = df_EC.reindex(timeIndex)#.ffill()

    print('df_EC loaded. Columns in this dataframe:')
    print(df_EC.keys())
    return df_EC

#--- 2. Load raw Storage flux data (no EC and without gapfilling)
def Read_LooStor(years,datapath):
    df_Stor         = []
    for iy in years:
        print('Loading %d'%iy)
        infilename  = '%s\Meteo\Storageflux\LooStor%02d.csv'%(datapath,np.mod(iy,100)) # just the storage data
        df_Stor       .append(pd.read_csv(infilename, index_col='YearDayTime', parse_dates=['YearDayTime'], date_parser=dateparse))
    df_Stor         =  pd.concat(df_Stor, axis=0, ignore_index=False)
    
    # Make sure the time index is complete, so all data frames have the same length.
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_Stor         = df_Stor.reindex(timeIndex)#.ffill()

    print('df_Stor loaded. Columns in this dataframe:')
    print(df_Stor.keys())
    return df_Stor

#--- 3. Load Combined data file (containing EC+Storage+Gapfilled)
def Read_LoodatGapfill(years,datapath):
    df_Comb     = []
    for iy in years:
        print('Loading %d'%iy)
        infilename            = '%s\Eddy\StorGapfilled\%4d\Loodat%02dGapfill.csv.out'%(datapath,iy,np.mod(iy,100))
       #df_Comb.append(pd.read_csv(infilename,  sep='\t',header=0,skiprows=[1,],date_parser=dateparse_Comb,parse_dates={'YearDayTime': ['Year','DoY','Hour']},index_col='YearDayTime'))
        df_Comb.append(pd.read_csv(infilename,  sep='\t',header=0,skiprows=[1,]))
    df_Comb                      = pd.concat(df_Comb, axis=0, ignore_index=False)
   
    df_Comb.index = pd.to_datetime((np.asarray(df_Comb['Year'], dtype='datetime64[Y]')-1970)+(np.asarray(         df_Comb['DoY' ]-1, dtype='timedelta64[D]')) +(np.asarray(np.floor(df_Comb['Hour']), dtype='timedelta64[h]')) + (np.asarray(60*np.mod(  df_Comb['Hour'],1), dtype='timedelta64[m]')))
    df_Comb.index.name = 'YearDayTime'

    # Make sure the time index is complete, so all data frames have the same length.
    df_Comb         = df_Comb[~df_Comb.index.duplicated()]
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_Comb         = df_Comb.reindex(timeIndex)
    
    #QA/QC
    df_Comb[(df_Comb > -10000) & (df_Comb < -9998)] = np.nan
    df_Comb['Ustar'][df_Comb['Ustar']==0.0000]  = np.nan
    df_Comb['NEE'  ][df_Comb['NEE'  ]==0.0000]  = np.nan
    df_Comb['LE'   ][df_Comb['LE'   ]==0.0000]  = np.nan
    df_Comb['H'    ][df_Comb['H'    ]==0.0000]  = np.nan
    print('df_Comb loaded. Columns in this dataframe:')
    print(df_Comb.keys())
    
    return df_Comb

#--- 4. # Loading Final datafiles (EC+Storage+Gapfilled, contains the same data as the combined datafiles, but only a selection.)	
def Read_Loobos_halfhourly(years,datapath):
    df_NEE          = []
    for iy in years:
        print('Loading %d'%iy)
        infilename  = '%s\Eddy\Loobos_halfhourly_%04d.csv'%(  datapath,iy)             # NEE as combined EC and Storage flux data, including gapfilling
        df_tmp      = pd.read_csv(infilename, index_col='DateTime'   , parse_dates=['DateTime'   ], date_parser=dateparse_Gapfilled)
      
        # The headers have changed over the years, make sure they are identical for all years.    
        df_tmp      = df_tmp.rename(columns={'Data.NEE_f'       : 'NEE_f'  })
        df_tmp      = df_tmp.rename(columns={'Data.H_f'         : 'H_f'    })
        df_tmp      = df_tmp.rename(columns={'Data.LE_f'        : 'LE_f'   })
        df_tmp      = df_tmp.rename(columns={'Data.Reco'        : 'Reco_f' })
        df_tmp      = df_tmp.rename(columns={'Data.GPP_f'       : 'GPP_f'  })
        df_tmp      = df_tmp.rename(columns={'NEE_f(umolm-2s-1)': 'NEE_f'  })
        df_tmp      = df_tmp.rename(columns={'H_f(Wm-2)'        : 'H_f'    })
        df_tmp      = df_tmp.rename(columns={'LE_f(Wm-2)'       : 'LE_f'   })
        df_tmp      = df_tmp.rename(columns={'NEE_fqc(-)'       : 'NEE_fqc'})
        df_tmp      = df_tmp.rename(columns={'H_fqc(-)'         : 'H_fqc'  })
        df_tmp      = df_tmp.rename(columns={'LE_fqc(-)'        : 'LEfqc'  })
        df_tmp      = df_tmp.rename(columns={'Reco(umolm-2s-1)' : 'Reco'   })
        df_tmp      = df_tmp.rename(columns={'GPP_f(umolm-2s-1)': 'GPP_f'  })
        df_tmp      = df_tmp.rename(columns={'Tair_f(degC)'     : 'Tair'   })
        df_tmp      = df_tmp.rename(columns={'Rg_f(Wm-2)'       : 'Rg_f'   })
        df_tmp      = df_tmp.rename(columns={'VPD(hPa)'         : 'VPD'    })
        df_tmp      = df_tmp.rename(columns={'Tsoil(degC)'      : 'Tsoil'  })
        df_tmp      = df_tmp.rename(columns={'rH(%)'            : 'rH'     })
        df_tmp      = df_tmp.rename(columns={'Ustar(m/s)'       : 'Ustar'  })
        df_tmp      = df_tmp.rename(columns={'R-ref(umolm-2s-1)': 'R-ref'  })
        df_tmp      = df_tmp.rename(columns={'E_0(degK)'        : 'E_0'    })
        df_NEE       .append(df_tmp)
        
    df_NEE          = pd.concat(df_NEE,  axis=0, ignore_index=False)
    if 'Unnamed: 0' in df_NEE.columns:
        df_NEE          = df_NEE.drop(columns='Unnamed: 0')

    # Make sure the time index is complete, so all data frames have the same length.
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_NEE          = df_NEE.reindex(timeIndex)#.ffill()

    print('Done')
    print('df_NEE loaded. Columns in this dataframe:')
    print(df_NEE.keys())
    return df_NEE

# 5. Loading Final Meteo files
def Read_Loobos_meteo(years,datapath):
    df_meteo        = []
    for iy in years:
        print('Loading %d'%iy)
        infilename  = '%s\Meteo\Final\Loodat%02dfinal.csv'%(datapath,np.mod(iy,100))
        df_meteo     .append(pd.read_csv(infilename, index_col='YearDayTime', parse_dates=['YearDayTime'], date_parser=dateparse))
    df_meteo        =  pd.concat(df_meteo,   axis=0, ignore_index=False)

    # Make sure the time index is complete, so all data frames have the same length.
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_meteo        = df_meteo.reindex(timeIndex)#.ffill()

    print('df_meteo loaded. Columns in this dataframe:')
    print(df_meteo.keys())	
    return df_meteo

# 6. Read Final Soil data files; ; soil moisture in % saturation, divide by porosity to get VWC in m3/m3
def Read_Loobos_soil(years,datapath):
    df_soil         = []
    for iy in years:
        infilename  = '%s\Soilmoist\Final\Loosoifinal%02d.csv'%(datapath,np.mod(iy,100))
        print('Loading %s...'%infilename)
        df_soil       .append(pd.read_csv(infilename, index_col='YearDayTime', parse_dates=['YearDayTime'], date_parser=dateparse))
    df_soil         =  pd.concat(df_soil,   axis=0, ignore_index=False)

    # remove duplicates, otherwise the reindex won't work
    df_soil         = df_soil[~df_soil.index.duplicated()]

    # Make sure the time index is complete, so all data frames have the same length.
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_soil         = df_soil.reindex(timeIndex)#.ffill()

    # Make sure all columns are in, regardless of the years selected.
    for key in ['ST-003','ST-020','ST-050','ST-100','Temp200','Temp201','Temp202','Temp203','Temp204','Temp205','Temp206','Temp207']:
        if not key in df_soil.keys():
            df_soil[key] = np.nan
    print('df_soil loaded. Columns in this dataframe:')
    print(df_soil.keys())
        
    #QA/QC
    I               = ((df_soil == 0.0) ) # -9999 is label for missing value
    df_soil[I]      = np.nan
    I = ( (df_soil.index >= np.datetime64('2005-04-11 14:30')) & (df_soil.index <= np.datetime64('2005-04-11 15:00')))
    df_soil.loc[I,['ST-003','ST-020','ST-050','ST-100']] = np.nan
    I = ( (df_soil.index >= np.datetime64('2014-11-29 20:00')) & (df_soil.index <= np.datetime64('2014-12-01 06:30')))
    df_soil.loc[I,['ST-003','ST-020','ST-050','ST-100']] = np.nan
    I = ( (df_soil.index >= np.datetime64('2014-12-02 00:00')) & (df_soil.index <= np.datetime64('2014-12-02 16:00')))
    df_soil.loc[I,['ST-003','ST-020','ST-050','ST-100']] = np.nan
    I = ( (df_soil.index >= np.datetime64('2017-12-04 00:00')) & (df_soil.index <= np.datetime64('2018-04-29 00:00')))
    df_soil.loc[I,['ST-003','ST-020','ST-050','ST-100']] = np.nan
    
    for key in ['Temp200','Temp201','Temp202','Temp203','Temp204','Temp205','Temp206','Temp207']:
        I = ( (df_soil[key].diff() < -5) | (df_soil[key].diff() > 5))
        df_soil.loc[I,key] = np.nan
        
    for key in ['ST-003','ST-020','ST-050','ST-100','Temp200','Temp201','Temp202','Temp203','Temp204','Temp205','Temp206','Temp207']:
        I = ( (df_soil[key] < -20) | (df_soil[key]> 25))
        df_soil.loc[I,key] = np.nan

    return df_soil
    
#--- 7. Load raw profile data (without storage and without gapfilling)
def Read_Loobos_profile(years,datapath):
    df_profile      = []
    for iy in years:
        print('Loading %d'%iy)
        infilename  = '%s\RawStorage\Loograd%02dfinal.csv'%(datapath,np.mod(iy,100)) # Profile data without storage or gapfilling
        df_profile    .append(pd.read_csv(infilename, index_col='YearDayTime', parse_dates=['YearDayTime'], date_parser=dateparse))
    df_profile      =  pd.concat(df_profile,   axis=0, ignore_index=False)

    # Make sure the time index is complete, so all data frames have the same length.
    timeIndex       = pd.DatetimeIndex(np.arange(np.datetime64('%4d-01-01 00:00:00'%years[0]),np.datetime64('%4d-12-31 23:59:59'%years[-1]),np.timedelta64(30,'m')).astype(np.datetime64))
    df_profile      = df_profile.reindex(timeIndex)#.ffill()

    print('df_profile loaded. Columns in this dataframe:')
    print(df_profile.keys())
    return df_profile