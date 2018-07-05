import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

result=[['122', '109', '2343', '220', '19'],
 ['15', '407', '37', '10', '102'],
 ['100', '100', '100', '100', '100'],
 ['113', '25', '19', '31', '112'],
 ['43', '219', '35', '33', '14'],
 ['132', '108', '256', '119', '14'],
 ['22', '48', '352', '51', '438']]

result = np.array(result, dtype=np.int)

fig=plt.figure(figsize=(5, 5), dpi=150)
ax1=fig.add_subplot(111, projection='3d')

xlabels = np.array(['10/11/2013', '10/12/2013', '10/13/2013',
                    '10/14/2013', '10/15/2013'])
xpos = np.arange(xlabels.shape[0])
ylabels = np.array(['A1','C1','G1','M1','M2','M3','P1'])
ypos = np.arange(ylabels.shape[0])

xposM, yposM = np.meshgrid(xpos, ypos, copy=False)

zpos=result
zpos = zpos.ravel()

dx=0.5
dy=0.5
dz=zpos

ax1.w_xaxis.set_ticks(xpos + dx/2.)
ax1.w_xaxis.set_ticklabels(xlabels)

ax1.w_yaxis.set_ticks(ypos + dy/2.)
ax1.w_yaxis.set_ticklabels(ylabels)

values = np.linspace(0.2, 1., xposM.ravel().shape[0])
colors = cm.rainbow(values)
ax1.bar3d(xposM.ravel(), yposM.ravel(), dz*0, dx, dy, dz, color=colors)
plt.show()