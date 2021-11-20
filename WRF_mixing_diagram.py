'''
19 November 2021
Ryann Ashley Wakefield 

This script computes mixing diagrams (MD) from NU-WRF output to compare control, dry and wet soil simulations
during the 2012 drought. The objective is to compare how boundary layer moisture and temperature evolved in response to 
changing soil moisture to determine if land-atmosphere feedbacks played a role in drought evolution. 
Soil moisture perturbation simulations were performed by coupling NASA's NU-WRF to its Land Information Systems
(LIS) to obtain coupled simulations. 
 
Mixing diagrams are used to quantify the total change in the moisture and energy budget of the daytime 
boundary layer, and to estimate the relative contributions of surface fluxes and entrainment
to that evolution. The methodology is introduced in Santanello et al. (2009)

A new methodology using mixed layer mean values (instead of 2-m) to represent PBL evolution is introduced in this script as part
of a new methodology in preparation by Wakefield et al. (2021)

*** reference: Joseph A. Santanello Jr., Christa D. Peters-Lidard, Sujay V. Kumar, Charles Alonge, and Wei-Kuo Tao, 2009: A Modeling and 
               Observational Framework for Diagnosing Local Land–Atmosphere Coupling on Diurnal Time Scales. 
               J. Hydrometeor, 10, 577–599. doi:10.1175/2009JHM1066.1
               
               Wakefield, R.A., J.A. Santanello, J.B. Basara: Sensitivity of daytime planetary boundary layer evolution to soil moisture 
               conditions during the 2012 flash drought. in preparation.    
'''



import numpy as np
from netCDF4 import Dataset
from wrf import getvar, ALL_TIMES, xy_to_ll, ll_to_xy
import xarray as xr
import sys
import os, re
import glob
import datetime
import metpy
import metpy.calc as mcalc
from metpy.units import units
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import rc
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D


#set up plotting properties
rc('font', weight='normal')
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['mathtext.default'] = 'regular'

plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.titleweight'] = 'normal'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 14

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.sans-serif'] = 'Computer Modern Roman'
plt.rcParams['font.weight'] = 'normal'	



'''
###
uses sub domain that is defined by max_dist_from_center and grid_box_spacing in execute_md_py.sh
the sub domain is for averaging profiles and eliminating noise if wrf output is high spatial resolution
sub domain is a square with length and width equal to 2*max_dist_from_center
###
'''

#define global constants
grav = units.Quantity(-9.8, 'm/(s**2)')
Lv = units.Quantity(2.5e6, 'J/kg')
Cp = units.Quantity(1004., 'J/kg/K')


#compute sub-domain indices
def get_gridbox_indices(igridx, n_km, grid_interval):
    x = lambda a, b, c:np.arange(a-b, a+b+c,c)
    return x(igridx, n_km, grid_interval)


'''
###
create function to extract wrf data using input gridbox locations and file list
###
'''

def get_wrfdata(wrflist, ix, iy):
    wrfvars = ['theta', 'pres', 'height_agl', 'TH2', 'PSFC', 'Q2', 'QVAPOR', 'LH', 'HFX', 'PBLH']
    
    #create list of meteorological variables averaged over the sub domain 
    datalist = [getvar(wrflist, x, timeidx=ALL_TIMES, method='cat').isel(south_north = iy, west_east = ix)\
            .mean(dim = ['south_north', 'west_east'], keep_attrs = True) for x in wrfvars]
    
    #create dictionary of wrf data arrays
    wrfdata = {x:y for (x,y) in zip(wrfvars, datalist)}

    #get time variables and compute difference between timesteps
    time = getvar(wrflist, 'times', timeidx=ALL_TIMES, method='cat')
    dt = ((time[1:]-time[:-1])/np.timedelta64(1,'s')).assign_attrs(units = 's') #ensure that timesteps are converted to seconds. 
    dt = dt.metpy.quantify().metpy.unit_array
    return wrfdata, time, dt

'''
###
define mixing diagram object so that mixing diagrams can be computed for several simulations
in this particular script, I want mixing diagrams for control, dry soil and wet soil simulations
so it is useful to have a mixing diagram object that corresponds to each
###
'''

class mixing_diagram_calculations:
    def __init__(self,wrfdata, ml_frac, sl_frac):
        #initialize arrays that correspond to mixing diagram curve
        self.lvq = np.zeros(shape = len(time))*np.nan
        self.cptheta = np.zeros(shape = len(time))*np.nan
        self.lvq2m = np.zeros(shape = len(time))*np.nan
        self.cptheta2m = np.zeros(shape = len(time))*np.nan
        self.pblmass = np.zeros(shape = len(time))*np.nan
        for i in range(len(time)):
            met_vars = ['pres', 'PBLH','height_agl','PSFC', 'theta', 'QVAPOR']
            p, pbl, h, psfc, theta, w = [wrfdata[x].isel(Time= i).metpy.quantify().metpy.unit_array for x in met_vars]
            p_layer = mcalc.get_layer(p, height = h, bottom = p[0], depth= pbl-h[0]) 
            pblp = p_layer[0][-1] ; 
            self.pblmass[i] = ((p_layer[0][-1]-psfc).to_base_units()/grav).magnitude

            '''
            define surface layer and mixed_layer depth
            we want to quantify the moisture and temperature of a mixed layer defined between
            sl_depth and ml_depth. Using ml_frac = 0.7 and sl_frac = 0.1, 
            If the PBL is 1000m deep, then sl_depth = 100m and ml_depth = 700 m. 
            this helps avoid superadiabatic surface layer and entrainment zone. 
            '''
            sl_depth = sl_frac*pbl ; ml_depth = (ml_frac*pbl)

            #compute bottom of layer for input into mixed layer mean function
            if sl_depth<h[0]:
                h_bot = h[0]
            else:
                h_bot = sl_depth
            
            #we want depth between lowest layer (above surface layer) and top layer (below entrainment zone)
            #layer depth is computed from bottom of our layer (surface layer) so to look between 0.1PBLH and 0.7PBLH
            # we need layer_depth = ml_depth-h_bot
            layer_depth = ml_depth-h_bot 
            thetalayer, wlayer = mcalc.mixed_layer(p, theta, w, height = h, bottom = h_bot, depth= layer_depth) 
            qlayer = mcalc.specific_humidity_from_mixing_ratio(wlayer).to('kg/kg') #convert to kg/kg just in case

            #mixed layer mixing diagram curve values. These represent the evolution of moisture and temperature
            #in the mixed convective boundary layer (CBL)
            self.cptheta[i] = (thetalayer*Cp).magnitude 
            self.lvq[i] = (qlayer*Lv).magnitude

            #2m mixing diagram curve values. Typically, the assumption is that these can represent
            #the evolution of the entire mixed CBL well, so we will extract them for comparison
            self.cptheta2m[i] = (Cp*(wrfdata['TH2'].isel(Time = i).metpy.quantify().metpy.unit_array)).magnitude
            sfcq = mcalc.specific_humidity_from_mixing_ratio(wrfdata['Q2'].isel(Time = i).metpy.quantify().metpy.unit_array)
            self.lvq2m[i] = (sfcq*Lv).magnitude

        #extract appropriate units for mixing diagram curves
        cpthetaunit = (thetalayer*Cp).units ; cptheta2munit = (Cp*wrfdata['TH2'].isel(Time = i).metpy.quantify().metpy.unit_array).units
        lvqunit = (qlayer*Lv).units ; lvq2munit = (sfcq*Lv).units
        massunit = ((p_layer[0][-1]-psfc).to_base_units()/grav).units

        #apply units to mixing diagram curve arrays and convert to kJ/kg
        self.cptheta = (self.cptheta*cpthetaunit).to('kJ/kg') ; self.cptheta2m = (self.cptheta2m*cptheta2munit).to('kJ/kg')
        self.lvq = (self.lvq*lvqunit).to('kJ/kg') ; self.lvq2m = (self.lvq2m*lvq2munit).to('kJ/kg')
        self.pblmass = self.pblmass*massunit

        #compute total change in moisture (lvq) and temperature(cptheta) from beginning to end of period
        #mixed layer
        self.dlvq = self.lvq[-1]-self.lvq[0]
        self.dcptheta = self.cptheta[-1]-self.cptheta[0]
        #2-m
        self.dlvq2m = self.lvq2m[-1]-self.lvq2m[0]
        self.dcptheta2m = self.cptheta2m[-1]-self.cptheta2m[0]


        '''
        change from t1 to t2 (i.e. 14 to 15 UTC or 15 to 16 UTC) 
        is related to the input of moisture and heat from surface fluxes
        therefore we take the average surface fluxes between t1 and t2 and normalize them to 
        the depth of the PBL (in mass per unit area) between t1 and t2
        '''

        #averages of surface flux and pbl mass
        lhbar = wrfdata['LH'].rolling(Time = 2).mean().dropna('Time').metpy.quantify().metpy.unit_array
        shbar = wrfdata['HFX'].rolling(Time = 2).mean().dropna('Time').metpy.quantify().metpy.unit_array
        pblmassbar = ((self.pblmass[1:]+self.pblmass[:-1])*0.5)


        #distribute fluxes throughout PBL to obtain values in units of KJ/kg
        self.sfx = (lhbar*dt/pblmassbar).to('kJ/kg')
        self.sfy = (shbar*dt/pblmassbar).to('kJ/kg')


		#derive entrainment estimate by subtracting cumulative surface fluxes from the
		#total daytime evolution
        self.entx = self.dlvq-np.nansum(self.sfx) #entrainment contribution to moisture budget
        self.enty = self.dcptheta-np.nansum(self.sfy) #entrainment contribution to heat budget
        
        #same as above, except using 2-m total evolution instead of mixed layer evolution
        self.entx2m = self.dlvq2m-np.nansum(self.sfx)#entrainment contribution to moisture budget
        self.enty2m = self.dcptheta2m-np.nansum(self.sfy)#entrainment contribution to heat budget
        
        #cumulative daytime surface flux contribution to budgets 
        self.sfxsum = np.nansum(self.sfx)#surface flux contribution to moisture budget
        self.sfysum = np.nansum(self.sfy)#surface fluxe contribution to heat budget
        


'''
###
import arguments for date, simulation hours, directory and location information. 
###
'''
print(sys.argv)
d0 = str(sys.argv[1])
yy, mm, dd = d0.split('-')
#make sure all month and day strings have leading zero
mm = mm.zfill(2)
dd = dd.zfill(2)
hr1 = sys.argv[2]
hr2 = sys.argv[3]
sim_lat = float(sys.argv[4])
sim_lon = float(sys.argv[5])
n_km = int(sys.argv[6]) #maximum distance in x and y directions from center point with units of number of grid boxes  
grid_int = int(sys.argv[7]) #grid_box_spacing in domain over which we're averaging. currently using with 1 km to correspond with simulation resolution
directory_path = sys.argv[8]

'''
###
create file lists for each simulation set 
###
'''
sim_typelist = ['ctrl', 'dry', 'wet']
flists = []
for i, sim_type in enumerate(sim_typelist):
    sim_type = sim_typelist[i].upper()
    directory = directory_path+yy+mm+dd+'/'+sim_type+'/'
    all_files = [os.path.basename(x) for x in sorted(glob.glob(directory+'*'))]    
    flists.append([directory+x for x in all_files\
                    if (int(re.findall('\d+', x)[-3]) >= int(hr1))& \
                    (int(re.findall('\d+', x)[-3]) <= int(hr2))])
    

'''
###
Obtain datasets for employing in wrf-python's get_var call so we can extract desired variables
###
'''
wrflist = [Dataset(x) for x in flists[0]]

'''
###
#find the nearest x,y indices for desired lat and lon
###
'''
x_y = ll_to_xy(wrflist[0], sim_lat, sim_lon) #returns corresponding x,y coordinate pairs, which are actually in the order y, x


'''
###
use x, y indices to define a sub domain over which we will average our meteorological data 
to compute the mixing diagrams. All wrf files in this analysis have the same size and domain, 
so we only need to calculate the iy and ix indices once. 
###
'''
iy = get_gridbox_indices(x_y[0], n_km, grid_int) 
ix = get_gridbox_indices(x_y[1], n_km, grid_int)


'''
###
Create mixing diagram objects for each simulation
###
'''
#get mixing diagram object for control simulation
wrfdata, time, dt = get_wrfdata(wrflist, ix, iy)
md1 = mixing_diagram_calculations(wrfdata, 0.7, 0.1)
del wrfdata
print(sim_typelist[0]+' computations done')


#get mixing diagram object for dry soil simulation
wrflist = [Dataset(x) for x in flists[1]]
wrfdata, time, dt = get_wrfdata(wrflist, ix, iy)
md2 = mixing_diagram_calculations(wrfdata, 0.7, 0.1)
del wrfdata
print(sim_typelist[1]+' computations done')

#get mixing diagram object for wet soil simulation
wrflist = [Dataset(x) for x in flists[2]]
wrfdata, time, dt = get_wrfdata(wrflist, ix, iy)
md3 = mixing_diagram_calculations(wrfdata, 0.7, 0.1)
del wrfdata
print(sim_typelist[2]+' computations done')



'''
###
###
PLOT THE MIXING DIAGRAMS
###
###
'''

###
####
####change this plotting path for your own directory. 
plotting_path = '...'
####
####
####




#
#list of md objects to be plotted
#
mdlist = [md1, md2,md3]
curvecolor = ['k', 'k', 'k'] #can make different colors for each simulation type if desired

#set axis limits. Make sure both have same "axint" to retain symmetry
xmin = 14
axint = 32
xmax = xmin+axint
ymin = 300
ymax = ymin+axint

#set up descriptive strings for plot 
titledate = datetime.datetime.strftime(datetime.datetime(int(yy), int(mm), int(dd)), '%d %b %Y')
savedate_string = datetime.datetime.strftime(datetime.datetime(int(yy), int(mm), int(dd)), '%d_%b_%Y')
location = str(sim_lat)+'N '+str(sim_lon)+'W'


#set up sizing/partitioning of space for subplots
gs = gridspec.GridSpec(28,28)
gslist = [gs[5:21, :8], gs[5:21, 10:18], gs[5:21, 20:28]]
titlelist = ['CTRL', 'DRY soil', 'WET soil']



fig = plt.figure(figsize = (18,8))
fig.suptitle(titledate+' '+hr1+'00 to '+hr2+'00 UTC\n'+location, fontsize = 22)
for i, md in enumerate(mdlist):
    ax = fig.add_subplot(gslist[i])
    ax.set_title(titlelist[i])
    vectorcolor = ['g', 'orange']
    txtlabel = ['surface flux','entrainment']
    
    
    #plot the mixed layer MD's
    ax.plot(md.lvq, md.cptheta, color = curvecolor[i], marker = 'o', linewidth = 2)
    #create mixing diagram vectors
    arrowbase = [(md.lvq[0], md.cptheta[0]),(md.lvq[0]+md.sfxsum, md.cptheta[0]+md.sfysum)]
    arrowtip = [(md.lvq[0]+md.sfxsum, md.cptheta[0]+md.sfysum), (md.lvq[-1], md.cptheta[-1])]
    nvectors = len(arrowtip)
    for j in range(0,nvectors):
        ax.annotate("", xy=arrowtip[j],xytext=arrowbase[j], arrowprops=dict(arrowstyle="-|>", mutation_scale=25,  color = vectorcolor[j], alpha = 1, linewidth = 3))
    
    
    #plot the 2-m MD's
    ax.plot(md.lvq2m, md.cptheta2m, color= curvecolor[i], marker = 'P', linewidth = 1)
    #create mixing diagram vectors
    arrowbase = [(md.lvq2m[0], md.cptheta2m[0]),(md.lvq2m[0]+md.sfxsum, md.cptheta2m[0]+md.sfysum)]
    arrowtip = [(md.lvq2m[0]+md.sfxsum, md.cptheta2m[0]+md.sfysum), (md.lvq2m[-1], md.cptheta2m[-1])]
    for j in range(0,nvectors):
        ax.annotate("", xy=arrowtip[j],xytext=arrowbase[j], arrowprops=dict(arrowstyle="-|>", mutation_scale=25,  color = vectorcolor[j], alpha = 1, linewidth = 1))

	#axis properties and scaling
    ax.xaxis.set_major_locator(mtick.MultipleLocator(5))
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(1))
    ax.yaxis.set_major_locator(mtick.MultipleLocator(5))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(1))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('$L_v q$'+' ['+r'$kJ\; kg^{-1}$'+']', fontsize = 20)
    ax.set_ylabel('$C_p$'+r'$\theta$'+' ['+r'$kJ\; kg^{-1}$'+']', fontsize = 20)
    ax.set_ylim(ymin, ymax)
    ax.xaxis.grid(which = 'both', linestyle = ':')
    ax.yaxis.grid(which = 'both', linestyle = ':')

    
line1 = Line2D([0], [0], label='mixed layer evolution', color='k', marker = 'o', linewidth=3, linestyle='-')
line2 = Line2D([0], [0], label='2-m evolution', color='k', marker = 'P', linewidth=1, linestyle='-')
line3 = Line2D([0], [0], label=txtlabel[0], color=vectorcolor[0], linewidth=4, linestyle='-')
line4 = Line2D([0], [0], label=txtlabel[1], color=vectorcolor[1], linewidth=4, linestyle='-')
handlelist = [line1, line2, line3, line4]
plt.legend(handles=handlelist, loc='center', ncol = 2, handlelength = 3, bbox_to_anchor=(-0.8, -0.5), fontsize = 18)
plt.savefig(plotting_path+savedate_string+'.'+hr1+'.to.'+hr2+'UTC.mixing_diagrams.jpeg', dpi = 300)
plt.show()