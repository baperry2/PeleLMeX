import cantera as ct
import numpy as np
import pandas as pd
from collections import OrderedDict
import ctable_tools as ctable

# Thermo Conditions
press = ct.one_atm
temp = 298.0
phi = 1.0
fuel = 'CH4'
oxid = 'O2:1.0, N2:3.76'
mechanism = 'drm19.cti'
trans = 'Mix'
progvars = ["CO2","H2O","CO","H2"]
ctable_specs = ["CO2","H2O","CO","H2","N2","O2","OH"]
outfile_prefix = 'prem_drm19_phi1_p1_t298'

# Flame Numerics
width = 0.1
loglevel = 1
ratio = 2
slope = 0.05
curve = 0.05
prune = 0.02
max_points = 10000

# Set up the flame
gas = ct.Solution(mechanism)
gas.set_equivalence_ratio(phi, fuel, oxid)
flame = ct.FreeFlame(gas, width=width)
flame.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
flame.set_max_grid_points(1,max_points)
flame.transport_model = trans

# Solve Flame
flame.solve()

# Get all desired data
data = pd.DataFrame()
data["X"] = flame.grid
data["T"] = flame.T
data["VEL"] = flame.velocity
data["RHO"] = flame.density_mass
data["DIFF"] = flame.thermal_conductivity / flame.cp_mass
data["VISC"] = flame.viscosity

specXdata = pd.DataFrame(flame.X.T,
                         columns=gas.species_names)
specYdata = pd.DataFrame(flame.Y.T,
                         columns=gas.species_names)

rxnrates = flame.net_production_rates.T * list(gas.molecular_weights)
specRRdata = pd.DataFrame(rxnrates,columns=gas.species_names)

# Compute prgress variable and its source
# ensure prog is monotonic
# source term 0 at min and max to ensure boundedness
data["PROG"] = specYdata[progvars].sum(axis=1)
prog = np.array(data["PROG"])
mono = np.min((prog[1:] - prog[:-1]) >= 0)
if not mono: print("WARNING: Nonmonotonic progress variable")
data["SRC_PROG"] = specRRdata[progvars].sum(axis=1)
data.loc[0,"SRC_PROG"] = 0.0
data.loc[len(flame.grid)-1,"SRC_PROG"] = 0.0

# Convert to CGS units
ctable.convert_chemtable_units(data)

# Make the chemtable
chemtable = data.drop(columns=["VEL","X","PROG"])
for spec in ctable_specs:
    chemtable["Y-"+spec] = specYdata[spec]
chemtable = pd.DataFrame(chemtable.values,
                         columns = chemtable.columns,
                         index = pd.MultiIndex.from_product([data["PROG"]],names=['PROG']))
ctable.write_chemtable_binary(outfile_prefix+'.ctb', chemtable, "1DFGM")
ctable.print_chemtable(chemtable)
#print(chemtable)

chemtable.drop(columns=["SRC_PROG"], inplace=True)
ctable.write_chemtable_binary(outfile_prefix+'_norxn.ctb', chemtable, "1DFGM")
ctable.print_chemtable(chemtable)

# Make the PMF dat files: Detailed Chem and Manifold
def write_dat_file(fname, df):
    with open(fname,'w') as fi:
        line1 = "".join(["VARIABLES ="]
                        + [' "{}"'.format(var.split(' ')[0]) for var in df.columns[:4]]
                        + [' "{}"'.format(var.upper()) for var in df.columns[4:]])
        fi.write(line1+'\n')
        line2 = " ZONE I={} FORMAT=POINT\n".format(df.shape[1])
        fi.write(line2)
        print('Reformated file has these variables:')
        print(line1)
        for idex, row in df.iterrows():
            linen = "".join(['{:<26.15g}'.format(x) for x in row])+'\n'
            fi.write(linen)

rename = {'VEL':'u', 'T':'temp','RHO':'rho'}
keepvars = ["X","T","VEL","RHO"]
df = data[keepvars].rename(columns=rename)
df2 = pd.DataFrame(data[["PROG","RHO"]].values, columns=['X0','XRHO'])

write_dat_file(outfile_prefix+'.dat', pd.concat([df,specXdata],axis=1))
write_dat_file(outfile_prefix+'_mani.dat', pd.concat([df,df2],axis=1))
