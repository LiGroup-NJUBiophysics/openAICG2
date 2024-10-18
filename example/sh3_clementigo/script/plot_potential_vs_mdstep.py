import numpy as np
import matplotlib.pyplot as plt

_kcal_to_kj = 4.1840
log = np.loadtxt('../output/sh3_clementigo.log',skiprows=10,usecols=(1,2),dtype=float)
mdstep = log[:,0]
potential = log[:,4] / _kcal_to_kj
fig,ax = plt.subplots(figsize=(4.6,4))
ax.plot(mdstep,potential,'-')
ax.set_xlabel('MD step')
ax.set_ylabel('Potential(kcal/mol)')
plt.savefig('../output/potential_vs_mdstep.png',dpi=300,bbox_inches='tight')
plt.show()
