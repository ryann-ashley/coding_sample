This script computes mixing diagrams (MD) from NU-WRF output to compare control, dry and wet soil simulations
during the 2012 drought. The objective is to compare how boundary layer moisture and temperature evolved in response to 
changing soil moisture to determine if land-atmosphere feedbacks played a role in drought evolution. 
Soil moisture perturbation simulations were performed by coupling NASA's NU-WRF to its Land Information Systems
(LIS) to obtain coupled simulations. This script reads wrf output and computes average meteorological variables over 
a user defined sub-domain. It was designed to read 1km horizontal resolution output, so it was necessary to perform
spatial averages to avoid too much noise and lack of representativeness when examining a single point to describe
PBL evolution over a larger region. 

Mixing diagrams are used to quantify the total change in the moisture and energy budget of the daytime 
boundary layer, and to estimate the relative contributions of surface fluxes and entrainment
to that evolution. The methodology is introduced in Santanello et al. (2009)

A new methodology using mixed layer mean values (instead of 2-m) to represent PBL evolution is introduced in this script as part
of a new methodology in preparation by Wakefield et al. (2021)

reference: Joseph A. Santanello Jr., Christa D. Peters-Lidard, Sujay V. Kumar, Charles Alonge, and Wei-Kuo Tao, 2009: A Modeling and 
           Observational Framework for Diagnosing Local Land–Atmosphere Coupling on Diurnal Time Scales. 
           J. Hydrometeor, 10, 577–599. doi:10.1175/2009JHM1066.1
               
           Wakefield, R.A., J.A. Santanello, J.B. Basara: Sensitivity of daytime planetary boundary layer evolution to soil moisture 
           conditions during the 2012 flash drought. in preparation.    

All files are held within a directory with format "wrfout_YYYYmmdd" and subdirectories "CTRL", "DRY" and "WET"
The file format itself is "wrfout_d01_YYYY-mm-dd_HH_00_00"

Execution of this script is performed using execute_md_py.sh where inputs for directory path, simulation hours, date, 
location and sub domain grid box spacing can be chosen. 
