#--- import packages

import numpy as np
import pandas as pd
from datetime import datetime, timedelta

#--- Model after optimisation

#Calculate canopy resistance and CO2 assimilation/respiration
def runAgs(df_profile_filter,df_Comb_filter,df_meteo_filter,df_EC_filter,fstr,gmin_input=0.25):
    
        #fluxes using the A-Gs scheme.
        co2_ppm   = df_profile_filter['CO2level1']
        epsi      = 1.   # epsilon
        sigma     = 5.67E-8
        Tair_K    = df_Comb_filter['Tair'] + 273.  # This is the air temperature
        Ts_K      = ((df_meteo_filter['L(o)corr'] / (epsi * sigma))**0.25) # This is the surface temperature, which should be used in the model
        Ts_C      = Ts_K - 273.
        conv_fac  = 101.3 / (8.314 * Tair_K)       # converstion factor, obtained via the ideal gas law. mol / m3
        co2_mgm3  = (co2_ppm * 44.01) * conv_fac   # concentration * conversion factor * molar mass CO2.  mgm3 = ppm * g/mol / mol/m3
        Ts        = Ts_C

        rho_1     = 1.225            # Density of air kg/m3

            # Fixed constants 
        Q10gm     = 2.0              # Parameter to calculate the mesophyll conductance
        Q10am     = 2.0              # Parameter to calculate max primary productivity
        Q10gamma  = 2                # Parameter to calculate the CO2 compensation concentration. (2 in IFS, 1.5 in DALES)

            # Reference temperatures calculation mesophyll conductance:
        T1gm      = 273 - 273        # Converted to degreesC
        T2gm      = 309 - 273        # IFS=309, DALES=301 (default)

            # Reference temperatues calculation max primary productivity:
        T1Am      = 273 - 273        # IFS=281, DALES=286 (C4))   Converted to degrees C
        T2Am      = 313 - 273        # IFS=309, DALES=301 

        nuco2q    = 1.6              # Ratio molecular viscosity water to carbon dioxide
        #gmin      = 0.25 / 1000.     # Cuticular (minimum) conductance. NOTE: = g_cu in IFS, with a factor 1000 difference (m/s)
        gmin      = gmin_input / 1000.     # Cuticular (minimum) conductance. NOTE: = g_cu in IFS, with a factor 1000 difference (m/s)
        ad        = 0.07             # Regression coefficient to calculate Cfrac (kpa-1)
        Kx        = 0.7              # Extinction coefficient PAR (mground / mleaf)

            # Maximum quantum use efficiency
        epsilon0  =  0.0144    # Maximum quantum use efficiency. mgCO2 / J PAR. Also named alpha

            # Vegetation specific constants
        gm298_umol    = 0.09                        # obtained from litature: Knauer et al. 2018: Effects of mesophyll conductance .....
        gm298         = gm298_umol / conv_fac       # converted to (mm/s)
        Ammax298      = 2.6                         # CO2 maximal primary productivity
        f0            = 0.89                        # Maximum value Cfrac 
        co2_comp298   = (42 * 44.01) * (1/24.45)    # from ppm to mg/m3. Got value 42 from the Atmospheric boundary layer book

        #LAI trees (m2 m-2)
        LAI           = 2.1                         # Obtained from data measurements in Loobos 2021.

            # Constant molar mass  
        constants_M_co2 = 44.01
        constants_M_air = 28.97

            # Calculate the CO2 compensation concentration (IFS eq. 8.92)
            # "The compensation   point Γ is defined as the CO2 concentration at which the net CO2 assimilation of a fully lit leaf becomes zero."

        co2_comp = co2_comp298 * Q10gamma ** ((Ts - 25) / 10) # equation 8.92. co2_comp = mg/m3.

            # Calculate the mesophyll conductance (IFS eq. 8.93)
            # "The mesophyll conductance gm describes the transport of CO2 from the substomatal cavities to the mesophyll cells where the carbon is fixed."

        gm       = (gm298 * Q10gm **((Ts -25)/10)) / ((1. + np.exp(0.3*(T1gm - Ts)))*(1. + np.exp(0.3*(Ts  - T2gm)))) 
        gm       = gm / 1000. # convert to m/s

            # Calculate CO2 concentration inside the leaf (Ci)
        fmin0    = gmin/nuco2q - (1./9.) * gm

            # Calculate the minimum value of Cfrac
        fmin     = gmin /(gmin +gm) # Formula from IFS
        # fmin    = -fmin0 + ((fmin0 **2) + 4* gmin/nuco2q *gm)**0.5 / (2. *gm) # formula from DALES

        VPD      = df_Comb_filter['VPD']/10     #Our measurement data from Loobos converted to kPa (/10). Ds in Dales

        VPDmax   = (f0 - fmin) /ad   # VPDmax in kPa. Dmax in Dalese

            # Calculate the fraction of the concentration inside the leaf in comparison with the surface of the leaf. 
        cfrac    = f0 * (1 - VPD/VPDmax) + fmin * (VPD/VPDmax) # f in IFS.
 
            # Absolute CO2 concentration (mg/m3)
        co2_abs  = co2_mgm3 

            # CO2 concentration in leaf (mg/m3)
        ci       = cfrac * (co2_abs - co2_comp) + co2_comp

            # Max gross primary production in high light conditions 
            #  line 439 / formula 8.94. Ammax is in mg/m2/s
        Ammax    = (Ammax298 * Q10am ** ((Ts - 25)/10)) / ((1. + np.exp(0.3*(T1Am - Ts)))*(1. + np.exp(0.3*(Ts - T2Am))))

            # Gross assimilation rate (Am, IFS eq. 8.97). In mg/m2/s
        Am       = Ammax * (1 - np.exp(-(gm *(ci - co2_comp) / Ammax))) 
 
            # Autotrophic dark respiration (IFS eq. 8.99). In mg/m2/s
        Rdark    = Am / 9
 
            # Photosynthetically active radiation (PAR), Ia
        PAR      = df_meteo_filter['PAR'] * 0.22 # measured Loobos data. Convert from umol m-2 s-1 to Jm-2s-1

            # Calculate e (maximum quantum use efficiency) Also named as alpha. mgCO2 / J PAR
        epsilon  = epsilon0 * (co2_abs - co2_comp)/(co2_abs + 2. * co2_comp) # Formula from DALES

            # calculate the gross primary productivity (mg/m2/s)            
        Ag       = (Am + Rdark) * (1 - np.exp((-epsilon * PAR)/(Am + Rdark))) - Rdark # Formula 8.98   #this appears to be an orphan variable -FB

            # Calculate upscaling from leaf to canopy: net flow CO2 into the plant (An) [-]   
        tempy    = epsilon * Kx * PAR / (Am + Rdark)

        def E1(x):
            # E1() approximation
                euler = 0.5772156649015329
                G     = np.exp(-euler)
                b     = (2*(1-G)/(G*(2-G)))**0.5
                h_inf = (1-G)*(G**2 - 6*G+12) / (3*G*(2-G)**2*b)
                q     = 20/47*x**(31/26.)**0.5
                h     = 1 / (1+x*x**0.5)
                E1    = np.exp(-x) / (G+(1-G)*np.exp(-x/(1-G))) * np.log(1+G/x-(1-G)/(h+b*x)**2)
                return E1

            # Calculate the net assimilation

                # 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant  
        E1_first    = E1(tempy * np.exp(-Kx*LAI))
        E1_second   = E1(tempy)
        An_canopy   = (Am + Rdark) * (1 - 1. / (Kx * LAI) * (E1_first - E1_second)) # code from DALES

                # 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level
        a1          = 1.0 / (1 - f0) #reminder: Maximum value Cfrac 
        Dstar       = VPDmax / (a1 * (f0 - fmin))

        #moving fstr to a input of the model
        #fstr        = 1.     # ranges from 0: values at wilting point, to 1: absence of moisture stress
        gcco2       =     LAI * (gmin / nuco2q + a1 * fstr * An_canopy / ((co2_abs - co2_comp) * (1. + VPD / Dstar))) # m/s

                # 3. calculate surface resistance for moisture and carbon dioxide
        #rs          = 1.67 / (gcco2) #wrong  
        rs          = 1./ (1.67 * gcco2) #correct++++
    
        rsCO2       = 1. / gcco2         # Surface resistance of CO2 in s/m

                # calculate the ra, aerodynamic resistance
        U           = df_EC_filter['Mea_Windsp']
        U_star      = df_EC_filter['U-star']
        ra          = U / (U_star**2)             # get the ra from the Loobos observations

        # 4.  calculate net flux of CO2 into the plant (An, mg/m2/s)
        An_final    = (co2_abs - ci) / (ra + rsCO2)   # should have as default a minus sign before the formula
        # The assimilation rate (A) is expressed as amount of CO2 assimilated per unit leaf area and time (mol m−2 s−1)
        # end of Jamie's code

        #I want to convert to umol carbon /m2/s . Molar weight of CO2 is 44.01g/mol
        An_umol = (An_final / 44.01 )*1000

        return(An_final, An_umol, rs, ra)

#version of Ags that outputs a dataframe instead of a collection of series
def runAgs_df(df_profile_filter,df_Comb_filter,df_meteo_filter,df_EC_filter,fstr,gmin_input=0.25):
    
        #fluxes using the A-Gs scheme.
        co2_ppm   = df_profile_filter['CO2level1']
        epsi      = 1.   # epsilon
        sigma     = 5.67E-8
        Tair_K    = df_Comb_filter['Tair'] + 273.  # This is the air temperature
        Ts_K      = ((df_meteo_filter['L(o)corr'] / (epsi * sigma))**0.25) # This is the surface temperature, which should be used in the model
        Ts_C      = Ts_K - 273.
        conv_fac  = 101.3 / (8.314 * Tair_K)       # converstion factor, obtained via the ideal gas law. mol / m3
        co2_mgm3  = (co2_ppm * 44.01) * conv_fac   # concentration * conversion factor * molar mass CO2.  mgm3 = ppm * g/mol / mol/m3
        Ts        = Ts_C

        rho_1     = 1.225            # Density of air kg/m3

            # Fixed constants 
        Q10gm     = 2.0              # Parameter to calculate the mesophyll conductance
        Q10am     = 2.0              # Parameter to calculate max primary productivity
        Q10gamma  = 2                # Parameter to calculate the CO2 compensation concentration. (2 in IFS, 1.5 in DALES)

            # Reference temperatures calculation mesophyll conductance:
        T1gm      = 273 - 273        # Converted to degreesC
        T2gm      = 309 - 273        # IFS=309, DALES=301 (default)

            # Reference temperatues calculation max primary productivity:
        T1Am      = 273 - 273        # IFS=281, DALES=286 (C4))   Converted to degrees C
        T2Am      = 313 - 273        # IFS=309, DALES=301 

        nuco2q    = 1.6              # Ratio molecular viscosity water to carbon dioxide
        #gmin      = 0.25 / 1000.     # Cuticular (minimum) conductance. NOTE: = g_cu in IFS, with a factor 1000 difference (m/s)
        gmin      = gmin_input / 1000.     # Cuticular (minimum) conductance. NOTE: = g_cu in IFS, with a factor 1000 difference (m/s)
        ad        = 0.07             # Regression coefficient to calculate Cfrac (kpa-1)
        Kx        = 0.7              # Extinction coefficient PAR (mground / mleaf)

            # Maximum quantum use efficiency
        epsilon0  =  0.0144    # Maximum quantum use efficiency. mgCO2 / J PAR. Also named alpha

            # Vegetation specific constants
        gm298_umol    = 0.09                        # obtained from litature: Knauer et al. 2018: Effects of mesophyll conductance .....
        gm298         = gm298_umol / conv_fac       # converted to (mm/s)
        Ammax298      = 2.6                         # CO2 maximal primary productivity
        f0            = 0.89                        # Maximum value Cfrac 
        co2_comp298   = (42 * 44.01) * (1/24.45)    # from ppm to mg/m3. Got value 42 from the Atmospheric boundary layer book

        #LAI trees (m2 m-2)
        LAI           = 2.1                         # Obtained from data measurements in Loobos 2021.

            # Constant molar mass  
        constants_M_co2 = 44.01
        constants_M_air = 28.97

            # Calculate the CO2 compensation concentration (IFS eq. 8.92)
            # "The compensation   point Γ is defined as the CO2 concentration at which the net CO2 assimilation of a fully lit leaf becomes zero."

        co2_comp = co2_comp298 * Q10gamma ** ((Ts - 25) / 10) # equation 8.92. co2_comp = mg/m3.

            # Calculate the mesophyll conductance (IFS eq. 8.93)
            # "The mesophyll conductance gm describes the transport of CO2 from the substomatal cavities to the mesophyll cells where the carbon is fixed."

        gm       = (gm298 * Q10gm **((Ts -25)/10)) / ((1. + np.exp(0.3*(T1gm - Ts)))*(1. + np.exp(0.3*(Ts  - T2gm)))) 
        gm       = gm / 1000. # convert to m/s

            # Calculate CO2 concentration inside the leaf (Ci)
        fmin0    = gmin/nuco2q - (1./9.) * gm

            # Calculate the minimum value of Cfrac
        fmin     = gmin /(gmin +gm) # Formula from IFS
        # fmin    = -fmin0 + ((fmin0 **2) + 4* gmin/nuco2q *gm)**0.5 / (2. *gm) # formula from DALES

        VPD      = df_Comb_filter['VPD']/10     #Our measurement data from Loobos converted to kPa (/10). Ds in Dales

        VPDmax   = (f0 - fmin) /ad   # VPDmax in kPa. Dmax in Dalese

            # Calculate the fraction of the concentration inside the leaf in comparison with the surface of the leaf. 
        cfrac    = f0 * (1 - VPD/VPDmax) + fmin * (VPD/VPDmax) # f in IFS.
 
            # Absolute CO2 concentration (mg/m3)
        co2_abs  = co2_mgm3 

            # CO2 concentration in leaf (mg/m3)
        ci       = cfrac * (co2_abs - co2_comp) + co2_comp

            # Max gross primary production in high light conditions 
            #  line 439 / formula 8.94. Ammax is in mg/m2/s
        Ammax    = (Ammax298 * Q10am ** ((Ts - 25)/10)) / ((1. + np.exp(0.3*(T1Am - Ts)))*(1. + np.exp(0.3*(Ts - T2Am))))

            # Gross assimilation rate (Am, IFS eq. 8.97). In mg/m2/s
        Am       = Ammax * (1 - np.exp(-(gm *(ci - co2_comp) / Ammax))) 
 
            # Autotrophic dark respiration (IFS eq. 8.99). In mg/m2/s
        Rdark    = Am / 9
 
            # Photosynthetically active radiation (PAR), Ia
        PAR      = df_meteo_filter['PAR'] * 0.22 # measured Loobos data. Convert from umol m-2 s-1 to Jm-2s-1

            # Calculate e (maximum quantum use efficiency) Also named as alpha. mgCO2 / J PAR
        epsilon  = epsilon0 * (co2_abs - co2_comp)/(co2_abs + 2. * co2_comp) # Formula from DALES

            # calculate the gross primary productivity (mg/m2/s)            
        Ag       = (Am + Rdark) * (1 - np.exp((-epsilon * PAR)/(Am + Rdark))) - Rdark # Formula 8.98   #this appears to be an orphan variable -FB

            # Calculate upscaling from leaf to canopy: net flow CO2 into the plant (An) [-]   
        tempy    = epsilon * Kx * PAR / (Am + Rdark)

        def E1(x):
            # E1() approximation
                euler = 0.5772156649015329
                G     = np.exp(-euler)
                b     = (2*(1-G)/(G*(2-G)))**0.5
                h_inf = (1-G)*(G**2 - 6*G+12) / (3*G*(2-G)**2*b)
                q     = 20/47*x**(31/26.)**0.5
                h     = 1 / (1+x*x**0.5)
                E1    = np.exp(-x) / (G+(1-G)*np.exp(-x/(1-G))) * np.log(1+G/x-(1-G)/(h+b*x)**2)
                return E1

            # Calculate the net assimilation

                # 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant  
        E1_first    = E1(tempy * np.exp(-Kx*LAI))
        E1_second   = E1(tempy)
        An_canopy   = (Am + Rdark) * (1 - 1. / (Kx * LAI) * (E1_first - E1_second)) # code from DALES

                # 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level
        a1          = 1.0 / (1 - f0) #reminder: Maximum value Cfrac 
        Dstar       = VPDmax / (a1 * (f0 - fmin))

        #moving fstr to a input of the model
        #fstr        = 1.     # ranges from 0: values at wilting point, to 1: absence of moisture stress
        gcco2       =     LAI * (gmin / nuco2q + a1 * fstr * An_canopy / ((co2_abs - co2_comp) * (1. + VPD / Dstar))) # m/s

                # 3. calculate surface resistance for moisture and carbon dioxide
        #rs          = 1.67 / (gcco2) #wrong  
        rs          = 1./ (1.67 * gcco2) #correct++++
    
        rsCO2       = 1. / gcco2         # Surface resistance of CO2 in s/m

                # calculate the ra, aerodynamic resistance
        U           = df_EC_filter['Mea_Windsp']
        U_star      = df_EC_filter['U-star']
        ra          = U / (U_star**2)             # get the ra from the Loobos observations

        # 4.  calculate net flux of CO2 into the plant (An, mg/m2/s)
        An_final    = (co2_abs - ci) / (ra + rsCO2)   # should have as default a minus sign before the formula
        # The assimilation rate (A) is expressed as amount of CO2 assimilated per unit leaf area and time (mol m−2 s−1)
        # end of Jamie's code

        #I want to convert to umol carbon /m2/s . Molar weight of CO2 is 44.01g/mol
        An_umol = (An_final / 44.01 )*1000

        result_df = pd.concat([An_final, An_umol,rs, ra], keys = ['an_mg/m2/s','an_umol/m2/s','rs','ra'],axis=1)
        return(result_df)

#Calculate canopy resistance and CO2 assimilation/respiration
def runAgs2(df_profile_filter,df_Comb_filter,df_meteo_filter,df_EC_filter,fstr,gmin_input=0.25):
    
        #fluxes using the A-Gs scheme.
        co2_ppm   = df_profile_filter['CO2level1']
        epsi      = 1.   # epsilon
        sigma     = 5.67E-8
        Tair_K    = df_Comb_filter['Tair'] + 273.  # This is the air temperature
        Ts_K      = ((df_meteo_filter['L(o)corr'] / (epsi * sigma))**0.25) # This is the surface temperature, which should be used in the model
        Ts_C      = Ts_K - 273.
        conv_fac  = 101.3 / (8.314 * Tair_K)       # converstion factor, obtained via the ideal gas law. mol / m3
        co2_mgm3  = (co2_ppm * 44.01) * conv_fac   # concentration * conversion factor * molar mass CO2.  mgm3 = ppm * g/mol / mol/m3
        Ts        = Ts_C

        rho_1     = 1.225            # Density of air kg/m3

            # Fixed constants 
        Q10gm     = 2.0              # Parameter to calculate the mesophyll conductance
        Q10am     = 2.0              # Parameter to calculate max primary productivity
        Q10gamma  = 2                # Parameter to calculate the CO2 compensation concentration. (2 in IFS, 1.5 in DALES)

            # Reference temperatures calculation mesophyll conductance:
        T1gm      = 273 - 273        # Converted to degreesC
        T2gm      = 309 - 273        # IFS=309, DALES=301 (default)

            # Reference temperatues calculation max primary productivity:
        T1Am      = 273 - 273        # IFS=281, DALES=286 (C4))   Converted to degrees C
        T2Am      = 313 - 273        # IFS=309, DALES=301 

        nuco2q    = 1.6              # Ratio molecular viscosity water to carbon dioxide
        #gmin      = 0.25 / 1000.     # Cuticular (minimum) conductance. NOTE: = g_cu in IFS, with a factor 1000 difference (m/s)
        gmin      = gmin_input / 1000.     # Cuticular (minimum) conductance. NOTE: = g_cu in IFS, with a factor 1000 difference (m/s)
        ad        = 0.07             # Regression coefficient to calculate Cfrac (kpa-1)
        Kx        = 0.7              # Extinction coefficient PAR (mground / mleaf)

            # Maximum quantum use efficiency
        epsilon0  =  0.0144    # Maximum quantum use efficiency. mgCO2 / J PAR. Also named alpha

            # Vegetation specific constants
        gm298_umol    = 0.09                        # obtained from litature: Knauer et al. 2018: Effects of mesophyll conductance .....
        gm298         = gm298_umol / conv_fac       # converted to (mm/s)
        Ammax298      = 2.6                         # CO2 maximal primary productivity
        f0            = 0.89                        # Maximum value Cfrac 
        co2_comp298   = (42 * 44.01) * (1/24.45)    # from ppm to mg/m3. Got value 42 from the Atmospheric boundary layer book

        #LAI trees (m2 m-2)
        LAI           = 2.1                         # Obtained from data measurements in Loobos 2021.

            # Constant molar mass  
        constants_M_co2 = 44.01
        constants_M_air = 28.97

            # Calculate the CO2 compensation concentration (IFS eq. 8.92)
            # "The compensation   point Γ is defined as the CO2 concentration at which the net CO2 assimilation of a fully lit leaf becomes zero."

        co2_comp = co2_comp298 * Q10gamma ** ((Ts - 25) / 10) # equation 8.92. co2_comp = mg/m3.

            # Calculate the mesophyll conductance (IFS eq. 8.93)
            # "The mesophyll conductance gm describes the transport of CO2 from the substomatal cavities to the mesophyll cells where the carbon is fixed."

        gm       = (gm298 * Q10gm **((Ts -25)/10)) / ((1. + np.exp(0.3*(T1gm - Ts)))*(1. + np.exp(0.3*(Ts  - T2gm)))) 
        gm       = gm / 1000. # convert to m/s

            # Calculate CO2 concentration inside the leaf (Ci)
        fmin0    = gmin/nuco2q - (1./9.) * gm

            # Calculate the minimum value of Cfrac
        fmin     = gmin /(gmin +gm) # Formula from IFS
        # fmin    = -fmin0 + ((fmin0 **2) + 4* gmin/nuco2q *gm)**0.5 / (2. *gm) # formula from DALES

        VPD      = df_Comb_filter['VPD']/10     #Our measurement data from Loobos converted to kPa (/10). Ds in Dales

        VPDmax   = (f0 - fmin) /ad   # VPDmax in kPa. Dmax in Dalese

            # Calculate the fraction of the concentration inside the leaf in comparison with the surface of the leaf. 
        cfrac    = f0 * (1 - VPD/VPDmax) + fmin * (VPD/VPDmax) # f in IFS.
 
            # Absolute CO2 concentration (mg/m3)
        co2_abs  = co2_mgm3 

            # CO2 concentration in leaf (mg/m3)
        ci       = cfrac * (co2_abs - co2_comp) + co2_comp

            # Max gross primary production in high light conditions 
            #  line 439 / formula 8.94. Ammax is in mg/m2/s
        Ammax    = (Ammax298 * Q10am ** ((Ts - 25)/10)) / ((1. + np.exp(0.3*(T1Am - Ts)))*(1. + np.exp(0.3*(Ts - T2Am))))

            # Gross assimilation rate (Am, IFS eq. 8.97). In mg/m2/s
        Am       = Ammax * (1 - np.exp(-(gm *(ci - co2_comp) / Ammax))) 
 
            # Autotrophic dark respiration (IFS eq. 8.99). In mg/m2/s
        Rdark    = Am / 9
 
            # Photosynthetically active radiation (PAR), Ia
        PAR      = df_meteo_filter['PAR'] * 0.22 # measured Loobos data. Convert from umol m-2 s-1 to Jm-2s-1

            # Calculate e (maximum quantum use efficiency) Also named as alpha. mgCO2 / J PAR
        epsilon  = epsilon0 * (co2_abs - co2_comp)/(co2_abs + 2. * co2_comp) # Formula from DALES

            # calculate the gross primary productivity (mg/m2/s)            
        Ag       = (Am + Rdark) * (1 - np.exp((-epsilon * PAR)/(Am + Rdark))) - Rdark # Formula 8.98   #this appears to be an orphan variable -FB

            # Calculate upscaling from leaf to canopy: net flow CO2 into the plant (An) [-]   
        tempy    = epsilon * Kx * PAR / (Am + Rdark)

        def E1(x):
            # E1() approximation
                euler = 0.5772156649015329
                G     = np.exp(-euler)
                b     = (2*(1-G)/(G*(2-G)))**0.5
                h_inf = (1-G)*(G**2 - 6*G+12) / (3*G*(2-G)**2*b)
                q     = 20/47*x**(31/26.)**0.5
                h     = 1 / (1+x*x**0.5)
                E1    = np.exp(-x) / (G+(1-G)*np.exp(-x/(1-G))) * np.log(1+G/x-(1-G)/(h+b*x)**2)
                return E1

            # Calculate the net assimilation

                # 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant  
        E1_first    = E1(tempy * np.exp(-Kx*LAI))
        E1_second   = E1(tempy)
        An_canopy   = (Am + Rdark) * (1 - 1. / (Kx * LAI) * (E1_first - E1_second)) # code from DALES

                # 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level
        a1          = 1.0 / (1 - f0) #reminder: Maximum value Cfrac 
        Dstar       = VPDmax / (a1 * (f0 - fmin))

        #moving fstr to a input of the model
        #fstr        = 1.     # ranges from 0: values at wilting point, to 1: absence of moisture stress
        gcco2       =     LAI * (gmin / nuco2q + a1 * fstr * An_canopy / ((co2_abs - co2_comp) * (1. + VPD / Dstar))) # m/s

                # 3. calculate surface resistance for moisture and carbon dioxide
        #rs          = 1.67 / (gcco2) #wrong  
        rs          = 1./ (1.67 * gcco2) #correct++++
    
        rsCO2       = 1. / gcco2         # Surface resistance of CO2 in s/m

                # calculate the ra, aerodynamic resistance
        U           = df_EC_filter['Mea_Windsp']
        U_star      = df_EC_filter['U-star']
        ra          = U / (U_star**2)             # get the ra from the Loobos observations

        # 4.  calculate net flux of CO2 into the plant (An, mg/m2/s)
        An_final    = (co2_abs - ci) / (ra + rsCO2)   # should have as default a minus sign before the formula
        # The assimilation rate (A) is expressed as amount of CO2 assimilated per unit leaf area and time (mol m−2 s−1)
        # end of Jamie's code

        #I want to convert to umol carbon /m2/s . Molar weight of CO2 is 44.01g/mol
        An_umol = (An_final / 44.01 )*1000

        return(An_final, An_umol, rs, ra, Ts_C)


def calc_LE(df_ET):
    #correcting outgoing Longwave: the sensor measures values between -20 and 10, but that it because the blackbody emission from the sensor itself (dependent on the temp of the sensor) is not taken into account.
    #thus we must take the output of the sensor and add the emitted longwave radiation of the sensor itself.
    #R_L(out)_corrected = R_L(out)_measured + R_L(out)_sensor, where R_L(out)_sensor = sigma*T(sensor)^4

    #constants:
    sigma = 5.67e-8 # W/m2/K4, Stefan-boltzmann constant
    epsilon = 1/0.98
    df_ET['L(o)_sensor'] = sigma*((df_ET['Te-L(o)']+273)**4)    #where Te-L(o) is in C
    df_ET['L(o)_corr'] = df_ET['L(o)'] + df_ET['L(o)_sensor'] # where L(o)_corr is corrected Longwave out (corrected for sensor's own temp)

    #Formula for leaf temp is: R_L(out)_corrected = epsilon * sigma * T_sfc^4 (where epsilon = 0.98-1.00, sigma = 5.67e-8 W/m2/K4, T_sfc in K)
    #rearrange formula to:
    df_ET['T_sfc'] = (df_ET['L(o)_corr'] / (epsilon*sigma)) ** (1/4)  # T_sfc output in K)
    df_ET['T_sfc_C'] = df_ET['T_sfc']-273

    #constants:
    e_sat_0 = 0.6107 # e_sat_0 = 0.6107 kPa or 610.7 Pa
    a = 7.5
    b = 237.3 # oC (geen typo)
    df_ET['e_sat'] = e_sat_0 * 10**(a*df_ET['T_sfc_C'] / (b+df_ET['T_sfc_C'])) #  T_sfc_C in oC

    #formal clausius-Clapeyron (aka August-Roche-Magnus) from wikipedia: e_sat = e_sat_0 * 10^( 17.6*Temp / 243+ Temp)  where e_sat is in hPa and Temp is in K

    #VPD(in Pa) = e_sat - e_act
    #VPD(in kg/kg) = q_sat - q_act
    Rd = 287 # J/kg K
    Rv = 462 # J/kg K
    # e = vapour pressure # in Pa of kPa
    # p = air pressure # in Pa of kPa

    #q = Rd/Rv * e/p
    #q_sat = Rd/Rv * e_sat/p
    df_ET['q_sat'] = Rd/Rv * df_ET['e_sat']/df_ET['p_kPa']

    #method 1 of calculating e_act: through VPD from dataset
    #note: this is giving negative values so I'm removing it for now.
    #VPD = e_sat - e_act -> e_act = e_sat - VPD
    df_ET['e_act_fromVPD'] = df_ET['e_sat'] - df_ET['VPD_adj']

    #q_act = Rd/Rv * e_act/p
    df_ET['q_act_fromVPD'] = Rd/Rv * df_ET['e_act_fromVPD']/df_ET['p_kPa'] #adding this to check

    #final step, subtract to get VPD for specific humidity
    #VPD_q = q_sat-q_act
    df_ET['VPDq_fromVPD']=df_ET['q_sat'] - df_ET['q_act_fromVPD'] #adding this to check if there's a substantial difference

    #method 2 of calculating e_act: through Rel Humidity from dataset
    # RH = e_act/e_sat *100 -> e_act = RH * e_sat /100
    df_ET['e_act_fromRH'] = (df_ET['rH']/100)*df_ET['e_sat']

    #q_act = Rd/Rv * e_act/p
    df_ET['q_act_fromRH'] = Rd/Rv * df_ET['e_act_fromRH']/df_ET['p_kPa']

    #final step, subtract to get VPD for specific humidity
    #VPD_q = q_sat-q_act
    df_ET['VPDq_fromRH']=df_ET['q_sat'] - df_ET['q_act_fromRH']

    #final step
    #ET = rho * Lv * VPD/rs
    #rho = 1.2 (approx value given by Michiel), Lv = 2260 kJ/kg (from google) Note: update to more accurate values when I can

    df_ET['ET'] = 1.2 * 2260000 * (df_ET['VPDq_fromRH']/(df_ET['rs']+df_ET['ra']))
    df_ET['ET_VPD'] = 1.2 * 2260000 * (df_ET['VPDq_fromVPD']/(df_ET['rs']+df_ET['ra']))

    return df_ET

#alternate calc_LE using L(o)corr from df_meteo instead of applying my own correction to L(o), because that's how run_Ags does it.
def calc_LE2(df_ET):
    #correcting outgoing Longwave: 
    
    #constants:
    sigma = 5.67e-8 # W/m2/K4, Stefan-boltzmann constant
    epsilon = 1/0.98
    #Formula for leaf temp is: R_L(out)_corrected = epsilon * sigma * T_sfc^4 (where epsilon = 0.98-1.00, sigma = 5.67e-8 W/m2/K4, T_sfc in K)
    #rearrange formula to:
    df_ET['T_sfc'] = (df_ET['L(o)corr'] / (epsilon*sigma)) ** (1/4)  # T_sfc output in K)
    df_ET['T_sfc_C'] = df_ET['T_sfc']-273

    #constants:
    e_sat_0 = 0.6107 # e_sat_0 = 0.6107 kPa or 610.7 Pa
    a = 7.5
    b = 237.3 # oC (geen typo)
    df_ET['e_sat'] = e_sat_0 * 10**(a*df_ET['T_sfc_C'] / (b+df_ET['T_sfc_C'])) #  T_sfc_C in oC

    #formal clausius-Clapeyron (aka August-Roche-Magnus) from wikipedia: e_sat = e_sat_0 * 10^( 17.6*Temp / 243+ Temp)  where e_sat is in hPa and Temp is in K

    #VPD(in Pa) = e_sat - e_act
    #VPD(in kg/kg) = q_sat - q_act
    Rd = 287 # J/kg K
    Rv = 462 # J/kg K
    # e = vapour pressure # in Pa of kPa
    # p = air pressure # in Pa of kPa

    #q = Rd/Rv * e/p
    #q_sat = Rd/Rv * e_sat/p
    df_ET['q_sat'] = Rd/Rv * df_ET['e_sat']/df_ET['p_kPa']

    #method 1 of calculating e_act: through VPD from dataset
    #note: this is giving negative values so I'm removing it for now.
    #VPD = e_sat - e_act -> e_act = e_sat - VPD
    df_ET['e_act_fromVPD'] = df_ET['e_sat'] - df_ET['VPD_adj']

    #q_act = Rd/Rv * e_act/p
    df_ET['q_act_fromVPD'] = Rd/Rv * df_ET['e_act_fromVPD']/df_ET['p_kPa'] #adding this to check

    #final step, subtract to get VPD for specific humidity
    #VPD_q = q_sat-q_act
    df_ET['VPDq_fromVPD']=df_ET['q_sat'] - df_ET['q_act_fromVPD'] #adding this to check if there's a substantial difference

    #method 2 of calculating e_act: through Rel Humidity from dataset
    # RH = e_act/e_sat *100 -> e_act = RH * e_sat /100
    df_ET['e_act_fromRH'] = (df_ET['rH']/100)*df_ET['e_sat']

    #q_act = Rd/Rv * e_act/p
    df_ET['q_act_fromRH'] = Rd/Rv * df_ET['e_act_fromRH']/df_ET['p_kPa']

    #final step, subtract to get VPD for specific humidity
    #VPD_q = q_sat-q_act
    df_ET['VPDq_fromRH']=df_ET['q_sat'] - df_ET['q_act_fromRH']

    #final step
    #ET = rho * Lv * VPD/rs
    #rho = 1.2 (approx value given by Michiel), Lv = 2260 kJ/kg (from google) Note: update to more accurate values when I can

    df_ET['ET'] = 1.2 * 2260000 * (df_ET['VPDq_fromRH']/(df_ET['rs']+df_ET['ra']))
    df_ET['ET_VPD'] = 1.2 * 2260000 * (df_ET['VPDq_fromVPD']/(df_ET['rs']+df_ET['ra']))

    return df_ET