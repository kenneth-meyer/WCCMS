import pandas as pd 
import math

path = "./circular_txt_data/"
out_path = "./result_circular/"
filename = "for_paraview.csv"
outfile = "for_paraview.csv"

df = pd.read_csv(path + filename)

length=len(df.index)
disp = []
# make magnitudes negative...
for i in range(0,length): # will this include the last column?
	disp.append(math.sqrt(df['u_x'][i]**2 + df['u_y'][i]**2 + df['u_z'][i]**2))
	if df['dot'][i] < 0:
		disp[i] = disp[i]*-1

df['displacement'] = disp

df.to_csv(out_path + outfile)