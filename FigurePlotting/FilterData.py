#Data filtering happens here (removing nonsensical data or out of bounds data, like -500C in temperature etc. Usually result from sensor error).

import numpy as np

def Filter_wrap(df_Comb,df_profile,df_meteo,df_EC,filterversion='default'):
    '''Function for filtering data based on excluding nonsensical and out-of-bounds values (sensor error).
    args are df_Comb,df_profile,df_meteo,df_EC,filterversion in that order.
    returns CO2,Locorr,VPD,Ustar,df_profile_filter,df_meteo_filter,df_Comb_filter,df_EC_filter in that order.'''
    
    # Make filter for GPP orginial data and not gapfilled
    #General filters
    I = ((df_Comb['GPP_fqc']==0)&(df_meteo['PAR']>0))
    #t = df_profile.index                                          
    #time = (t < np.datetime64('2013-05-08')) | (t > np.datetime64('2013-06-01'))
    
    # Filter for CO2 data
    CO2 = (df_profile['CO2level1'] > 300)
    
    # Filter for L(o)corr data
    Locorr= (df_meteo['L(o)corr']>0) 
    
    # Filter for VPD data
    VPD = (df_Comb['VPD']>=0)
    
    # Filter for U-star
    Ustar = (df_EC['U-star']>=0)
    
    # Combine all filters
    filter = I & CO2 & Locorr & VPD & Ustar
    
    #Column 'CO2' is input from df_profile
    #df_profile_CO2 = df_profile[CO2]
    #df_profile_filter = df_profile_CO2[I]
    df_profile_filter = df_profile[filter]
    
    #Column 'L(o)corr' and 'PAR' are inputs from df_meteo,
    #df_meteo_CO2 = df_meteo[CO2]
    #df_meteo_filter = df_meteo_CO2[I]
    df_meteo_filter = df_meteo[filter]
    
    #Columns 'VPD' and 'Tair' are inputs from df_Comb
    #df_Comb_CO2 = df_Comb[CO2]
    #df_Comb_filter = df_Comb_CO2[I]
    df_Comb_filter = df_Comb[filter]
    
    # Columns 'Mea_Windsp' and 'U-star' are inputs from df_EC
    #df_EC_CO2 = df_EC[CO2]
    #df_EC_filter = df_EC_CO2[I]
    df_EC_filter = df_EC[filter]
    
    return CO2,Locorr,VPD,Ustar,df_profile_filter,df_meteo_filter,df_Comb_filter,df_EC_filter

def Filter_GPP_LE_NEE_VPD(df_Comb,fqc=None): #a special version of the previous function that only takes input df_Comb, for when you're only interested in GPP,LE,NEE,VPD
    '''Function for filtering data based on excluding nonsensical and out-of-bounds values (sensor error).
    args are df_Comb,fqc (quality flag for GPP).
    if fqc is not None, it will apply  (df_Comb['GPP_fqc']==fqc) over the entire dataset (removing LE, NEE, VPD values as well). 
    Use fqc = 0 to select original data only, not gapfiled'''

    #this function filters data a little bit more intelligently that the one above because it just sets the values to NaN instead of cross-eliminating any row where ANY row has bad values
    import numpy as np
    
    # Filter for VPD data
    VPD = df_Comb['VPD']
    df_Comb['VPD']=np.where(VPD>=0,VPD,np.NaN)

    # Filter for GPP data
    GPP = df_Comb['GPP_f']
    df_Comb['GPP_f']=np.where(GPP>-20,GPP,np.NaN)
    GPP = df_Comb['GPP_f']
    df_Comb['GPP_f']=np.where(GPP<100,GPP,np.NaN)

    
    # Filter for LE data
    LE = df_Comb['LE']
    df_Comb['LE']=np.where(LE>-300,LE,np.NaN)

    # Filter for NEE data
    NEE = df_Comb['NEE']
    df_Comb['NEE']=np.where(NEE>-20,NEE,np.NaN)
    
    if fqc!=None:
        # optional: Make filter for GPP orginial data and not gapfilled
        #General filters
        I = (df_Comb['GPP_fqc']==fqc)  #set to 0 for original data only, not gapfilled
        df_Comb = df_Comb[I]

    return df_Comb
    