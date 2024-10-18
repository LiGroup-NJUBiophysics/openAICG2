import sys
#sys.path.append('../../../')
from openaicg2.forcefield.functionterms.dihedral_terms import aicg_dihd_term

import openmm as mm
import numpy as np
import pandas as pd
from openmm import unit
from openmm import app
import mdtraj as md
import matplotlib.pyplot as plt

# get energy and force for selected potential energy function
def get_energy_and_force(simulation,positions,groups_index):
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy=True,groups={groups_index})
    force = simulation.context.getState(getForces=True,groups={groups_index}).getForces(asNumpy=True)
    #force = [np.sqrt(f.dot(f)) for f in force]
    return state.getPotentialEnergy()._value,force
def aicg_function(dihd,para):
    return -para[0] * np.exp(-(dihd - para[1])**2/(2 * para[2]**2))

# Create a simple system
box_len = 5*unit.nanometer
N_particle = 4
T = 300
system = mm.System()
system.setDefaultPeriodicBoxVectors([box_len,0,0],[0,box_len,0],[0,0,box_len])
for _ in range(N_particle):
    system.addParticle(137*unit.amu)

ang = 177
x = 0 + np.cos(ang*np.pi/180)
z = 0 + np.sin(ang*np.pi/180)
positions = np.array([[0.38,0,0],[0,0,0],[0,0.38,0],[x,0.38,z]])
positions = positions * unit.nanometer

# make topology 
num_residue = 1
num_atom_per_resi = 4
top = app.Topology()
top.addChain('CCA')
for i in range(num_residue):
    backbone = top._chains[-1]
    top.addResidue('CGA',backbone)
    residue = list(top.residues())[-1]
    atom1 = top.addAtom('CA',app.Element.getBySymbol('C'),residue)
    for j in range(num_atom_per_resi-1):
        atom2 = top.addAtom('CA',app.Element.getBySymbol('C'),residue)
        top.addBond(atom1,atom2)
        atom1 = atom2

# parameter 
k_gaussian = 1.23
sigma_gaussion = 0.15
nat_theta = 170*np.pi/180
_A_to_nm_ = 0.1
_kcal_to_kj_ = 4.1840

# make a dihedral parameter table
num_dihd = 1
pd_aicg_dihd_idx_para_all = []
for i in range(num_dihd):
    i_dihd_idx = [[0, 1, 2, 3]]
    dihd_para = [[k_gaussian*_kcal_to_kj_,nat_theta,sigma_gaussion]]
    pd_dihd_idx = pd.DataFrame(i_dihd_idx,columns=['a1','a2','a3','a4'])
    pd_dihd_para = pd.DataFrame(dihd_para,columns=['epsilon','natdihd','width'])
    pd_dihd_idx_para = pd.concat([pd_dihd_idx,pd_dihd_para],axis=1)
    pd_aicg_dihd_idx_para_all.append(pd_dihd_idx_para)
pd_aicg_dihd_idx_para_all = pd.concat(pd_aicg_dihd_idx_para_all,axis=0)

aicg_dihd_force_p = aicg_dihd_term(pd_aicg_dihd_idx_para_all,force_group=1)
system.addForce(aicg_dihd_force_p)
integrator = mm.LangevinIntegrator(T*unit.kelvin,1.0/unit.picosecond,2.0*unit.femtosecond)
# create a simulatoin 
simulation = app.Simulation(top,system,integrator)

dihd_ang = np.linspace(-np.pi,np.pi,100)
energy_p = []
dihd_all = []
for i in range(len(dihd_ang)):
    ang = -dihd_ang[i]
    x = 0 + np.cos(ang)
    z = 0 + np.sin(ang)
    positions = np.array([[0.38,0,0],[0,0,0],[0,0.38,0],[x,0.38,z]])
    positions = positions * unit.nanometer
    ener_periodic, force_p = get_energy_and_force(simulation,positions,1)
    traj = md.Trajectory(positions,top)
    dihd = md.compute_dihedrals(traj,[[0, 1, 2, 3]])
    energy_p.append(ener_periodic)
    dihd_all.append(dihd[0][0])

dihd_all = np.array(dihd_all)
fig,ax = plt.subplots(figsize=(4.6,4))
ax.plot(dihd_all*180/np.pi,np.array(energy_p)/_kcal_to_kj_,'o',markersize=5,color='r',label='periodic deal',alpha=0.3)
ax.set_xlabel(r'$\theta$',fontsize=20)
ax.set_ylabel('Energy (kcal/mol)',fontsize=20)
plt.savefig('result/aicg_dihd_term.png',dpi=600,bbox_inches='tight')
plt.show()
