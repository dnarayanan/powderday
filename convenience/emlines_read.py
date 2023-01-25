import numpy as np

#define the emission lines filename
fname = "/blue/narayanan/prerakgarg/pd_runs/m25n512/rinner_fix/snap305_pdva_miles/emlines.galaxy99.txt"
data = np.genfromtxt(fname,skip_header=1)
data = np.atleast_2d(data)
data_wav = np.genfromtxt(fname,max_rows=1)

#identify the line that you want to isolate -- here, O3/5007
O3_id = np.abs(5008.24 - data_wav).argmin()

O3_lum = 0
for j in range(len(data)):
    O3_lum += data[j][O3_id]

