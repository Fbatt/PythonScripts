# Settings
Username   = 'molen050'
years      = range(2012,2014)#(1997,2021)

#--- import packages
import os
datapath   = os.path.join('C:\\','Users',Username, 'Onedrive - Wageningen University & Research','Documents - LoobosTeamsite','DataShare')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import cufflinks as cf
cf.go_offline()
cf.set_config_file(offline=False, world_readable=True)
from datetime import datetime, timedelta
import sys
sys.path.insert(0, os.path.join(datapath,'PythonScripts'))
from Loobos_Toolbox import Read_LoobosEddFinal, Read_LooStor, Read_LoodatGapfill, Read_Loobos_halfhourly, Read_Loobos_meteo, Read_Loobos_soil, Read_Loobos_profile



if not 'progress' in globals(): progress = list()
if not 'dataloaded' in progress:
    print('datapath is set to %s'%datapath)

    # Read files
    df_EC           = Read_LoobosEddFinal    (years,datapath)
#   df_Stor         = Read_LooStor           (years,datapath)
#   df_Comb         = Read_LoodatGapfill     (years,datapath)
#   df_NEE          = Read_Loobos_halfhourly (years,datapath)
#   df_meteo        = Read_Loobos_meteo      (years,datapath)
#   df_soil         = Read_Loobos_soil       (years,datapath)
#   df_profile      = Read_Loobos_profile    (years,datapath)
#   progress.append('dataloaded')

#df_EC loaded. Columns in this dataframe:
#Index(['Doy', 'Dtime', 'Flx_Tsonic', 'Flx_Lo-H2O', 'Flx_Lo-CO2', 'Qf_Tsonic',
#       'Qf_Lo-H2O', 'Qf_Lo-CO2', 'Mea_Windsp', 'Mea_Tsonic', 'Mea_Lo-H2O',
#       'Mea_Lo-CO2', 'U-star', 'Z-over-L', 'Wind-Dir', '80PercFlux'],
#      dtype='object')


#--- Plot CO2 fluxes from different files on top of each other
df_plot         = pd.concat([df_EC['Flx_Lo-CO2'],df_Stor['TotalCO2flux'],df_Comb['NEE_orig']],axis=1,sort=False)
df_plot.columns = ['EC','Storage','Combined']
f1              = df_plot.iplot(asFigure=True, layout=dict(yaxis=dict(title='NEE (umol/m2/s) '), xaxis=dict(title='time')), width=2)
f1.show()

# plot soil temperature
df_plot         = df_soil[['ST-003','ST-020','ST-050','ST-100']]
f2              = df_plot.iplot(asFigure=True,layout=dict(yaxis=dict(title='Tsoil (oC) '), xaxis=dict(title='time')), width=2)
f2.show()

# plot soil moisture
porosity = 0.35
df_plot = df_soil[['SM-003','SM-020','SM-050','SM-100']]/porosity
f3 = df_plot.iplot(asFigure=True,layout=dict(yaxis=dict(title='SMC (m3/m3) '), xaxis=dict(title='time')), width=2)
f3.show()

# plot NEE as a function of Hour
#df_Comb.plot(x='Hour',y='NEE',style='o')
#plt.show()

# plot mean diurnal cycle of NEE
f5 = df_Comb.groupby(df_Comb['Hour'])['NEE'].mean().iplot(asFigure=True,layout=dict(yaxis=dict(title='NEE (umol/m2/s) '), xaxis=dict(title='time')), width=2)
f5.show()

# Convert pandas dataframe to numpy array
np_NEE = np.asarray(df_NEE)
t = np.asarray(df_NEE.index)

# Plot data as numpy array using matplotlib.pyplot (as plt)
NEE = np_NEE[:,0]  # Column order remains the same.
NEE[(NEE < -100) | (NEE > 50)] = np.nan # remove outliers

f = plt.figure(1)
ax = f.add_subplot(111)
ax.plot(t,NEE) 
ax.set_xlabel('time')
ax.set_ylabel('NEE (umol/m2/s)')


