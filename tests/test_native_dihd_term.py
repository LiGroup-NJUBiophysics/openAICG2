import numpy as np
import matplotlib.pyplot as plt
def native_dihd_func(dihd,para):
    return (para[0] * (1 - np.cos(dihd-para[1])) + para[2] * (1 - np.cos(3 * (dihd-para[1])))).astype('float64')

def native_dihd_func_p(dihd,para):
    d_dihd = dihd - para[1]
    ddihd_periodic = d_dihd - np.floor((d_dihd+np.pi)/(2*np.pi))*(2*np.pi)
    return (para[0] * (1 - np.cos(ddihd_periodic)) + para[2] * (1 - np.cos(3 * ddihd_periodic))).astype('float64')

theta = np.linspace(-np.pi,np.pi,500)

value = native_dihd_func(theta,[1,-170,0.5])
valuep = native_dihd_func_p(theta,[1,-170,0.5])
fig,ax = plt.subplots(figsize=(4.6,4))
ax.plot(theta*180/np.pi,value,'*',color='k',label='no periodic')
ax.plot(theta*180/np.pi,valuep,'--',color='r',label='with periodic')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel('Value')
plt.legend()
plt.show()
