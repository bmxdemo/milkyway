from vlsr import vlsr
import bmxobs as bmx
import numpy as np 
import matplotlib.pyplot as plt

data = bmx.BMXObs("/home/chandrahas/Desktop/BMX/pas/210903_0000", channels="111")

vel = vlsr(data.ra*180/np.pi, 40.8, data.mjd)
v_rad = vel[0]+vel[1]+vel[2]
print(v_rad)
print(np.isnan(np.sum(v_rad)))
plt.plot(np.arange(len(v_rad)),v_rad)
plt.show()