import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
from scipy.stats import binned_statistic 

parser = argparse.ArgumentParser(
    prog = 'Statistics Plotting Utility',
    description = 'Reads PeleLMeX 2D plotfiles and plots conditional statistics',
    epilog = 'A work in rogress.')

# Default values for parsed arguments
filepattern = "stats/zoverD010/plt029*"
pltvars = ['temp','z_velocity','ZMIX','PROG','MANI_Y-CO']
condvars = ['RADIUS','ZMIX']
condvar_vals= {'RADIUS':np.linspace(0,0.0225,101),
               'ZMIX':np.linspace(0,0.5,78)}

def get_exp_data(varname, data, dattype='mean'):
    translation = {'RADIUS':'xmm', 'ZMIX':'Fblgr', 'temp':'Tray', 'MANI_Y-CO':'COLIF'}
    conversion = {'RADIUS':1000, 'ZMIX':(1.0/4.6)}
    if varname in translation.keys():
        var = translation[varname]
    else:
        var = varname
    if varname in conversion.keys():
        factor = conversion[varname]
    else:
        factor = 1
    return np.abs(data[var +' '+ dattype]) / factor
    
# Parse arguments
parser = argparse.ArgumentParser(
    prog = 'Statistics Plotting Utility',
    description = 'Reads PeleLMeX 2D plotfiles and plots conditional statistics',
    epilog = 'A work in rogress.')
parser.add_argument("-f","--filepattern",help="will read all files that match this pattern",default=filepattern)
parser.add_argument("-d","--depvars",help="dependent variables", default=pltvars, nargs='*')
parser.add_argument("-i","--indvars",help="independent variables", default=condvars, nargs='*')
parser.add_argument("-s","--save",help="flag to save plots", action='store_true')
args = parser.parse_args()

expdata = {}
datafile = '/projects/hpacf/bperry/sydney/data5GP/Reacting_Species_Mean_RMS-5GP/Radial_Favre/FJ200-5GP-Lr300-59-xDDISTANCE0_Radial_Favre.csv'
idx = args.filepattern.index('zoverD')+len('zoverD') +1
expdata['RADIUS'] = pd.read_csv(datafile.replace('DISTANCE',args.filepattern[idx:idx+2]), dtype=np.double)
datafile = '/projects/hpacf/bperry/sydney/data5GP/Reacting_Species_Mean_RMS-5GP/Mixture_fraction_Favre/FJ200-5GP-Lr300-59-xDDISTANCE0_Fblgr(+x)_Favre.csv'
expdata['ZMIX'] = pd.read_csv(datafile.replace('DISTANCE',args.filepattern[idx:idx+2]), dtype=np.double)
print(expdata['RADIUS'])
print(expdata['ZMIX'])


# Load the simulation plot files
ds = yt.load(args.filepattern)
allvars = list(set(args.depvars + args.indvars))

# Setup function to compute derived variables
def add_function(field, data):
    name = field.name[1]
    if name == 'ZMIX':
        return data['Y(ZMIX)']
    elif name == 'PROG':
        return data['Y(PROG)']
    elif name == 'RADIUS':
        vals = np.sqrt(data.fcoords[:,0]**2 + data.fcoords[:,1]**2)
        #vals.units = 'dimensionless' 
        return vals
    else:
        print(name)
        raise RuntimeError("Invalid variable requested: " + name)

# Read from the data files
data_accum = []
for dd in ds:
    # Try to derive any field that is not already present
    fields = [field[1] for field in dd.field_list]
    for field in allvars:
        if field not in fields:
            print(field)
            dd.add_field(('gas', field), add_function, 'cell')

    # Extract the data we need. Take from finest level, covered cells are ignored
    dt = dd.all_data()
    dt.min_level = dd.max_level
    #dt.max_level = 4
    npoints = len(dt[allvars[0]])
    print("Reading {} points from file: {}".format(npoints,dd) )
    data_accum.append(np.stack([np.array(dt[var]) for var in allvars]))

alldata = pd.DataFrame( np.concatenate(data_accum, axis=1).T, columns=allvars)

print(alldata)

# Take staistics
for condvar in args.indvars:
    if condvar in condvar_vals.keys():
        bins = condvar_vals[condvar]
    else:
        bins = 100

    means, bins, _ = binned_statistic(alldata[condvar],
                             [alldata[var] for var in args.depvars],
                             bins=bins, statistic='mean')
    stds, bins, _ =  binned_statistic(alldata[condvar],
                            [alldata[var] for var in args.depvars],
                             bins=bins, statistic='std')
    bincents = 0.5*(bins[1:] + bins[:-1])

    for ii,var in enumerate(args.depvars):
        txt = var + '_vs_' + condvar
        plt.figure(txt)
        plt.clf()
        plt.plot(bincents, means[ii], 'k-')
        plt.plot(bincents, stds[ii], 'r-')
        if ('velocity' not in var) and (var != 'PROG'):
            plt.plot(get_exp_data(condvar, expdata[condvar]), get_exp_data(var, expdata[condvar]), 'ko' )
            plt.plot(get_exp_data(condvar, expdata[condvar]), get_exp_data(var, expdata[condvar], 'RMS'), 'ro' )
        plt.xlabel(condvar)
        plt.ylabel(var)
        if args.save:
            path, pattern = os.path.split(args.filepattern)
            path = os.path.join('plots',path)
            if not os.path.exists(path): os.makedirs(path)
            figname = os.path.join(path, pattern[:6] + '_' + txt+'.pdf')
            print('Saving: ' + figname)
            plt.savefig(figname)

if not (args.save) :
    plt.show()
print(alldata)
