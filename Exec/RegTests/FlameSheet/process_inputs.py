import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cmlm.ctable_tools as ctt
from scipy.interpolate import NearestNDInterpolator

data = pd.read_csv('Ascent_Rig_allCCW_newDome_mfi_po_react_sbes_2019r2_5pre_3FTTsamp_prmx_exit.ascii')
data.columns = [col.strip() for col in data.columns]
# variables that will be lookup table inputs
invars = ['y-coordinate', 'z-coordinate']
# variables that will be lookup table ouputs
outvars =['mean-f','mean-x-velocity','mean-y-velocity','mean-z-velocity']

# grid definition for lookup table
ngrid = 801
yvals = np.linspace(-0.025, 0.025, ngrid);
zvals = np.linspace(-0.025, 0.025, ngrid);
gridy, gridz = np.meshgrid(yvals, zvals)

# interpolate unstructured grid onto structured table (nearest neighbor for now)
interper = NearestNDInterpolator(data[invars], data[outvars])
output = interper(gridy.flatten(), gridz.flatten())

# reverse z and y coordinate order -> order gets flipped when read in Pele
df_index = pd.MultiIndex.from_arrays([gridy.flatten(), gridz.flatten()], names=invars)
ctb_data = pd.DataFrame(index=df_index, columns=outvars, data=output)
ctt.write_chemtable_binary("bc_inflow_data.ctb", ctb_data, "TabulatedData", "Pele")

plt.figure()
plt.scatter(data['y-coordinate'], data['z-coordinate'], c=data['mean-f'], s=1)
plt.figure()
plt.contourf(np.array(ctb_data.index.get_level_values('y-coordinate')).reshape((ngrid, ngrid)),
             np.array(ctb_data.index.get_level_values('z-coordinate')).reshape((ngrid, ngrid)),
             np.array(ctb_data['mean-f']).reshape((ngrid, ngrid)),
             levels = 25)
plt.show()
