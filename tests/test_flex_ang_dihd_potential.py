import sys
from openaicg2.forcefield.functionterms.angle_terms import flex_angle_term
from openaicg2.forcefield.functionterms.dihedral_terms import flex_dihd_term
from openaicg2.utils import parser_flp_para
import openmm as mm
import numpy as np
import pandas as pd
from openmm import unit
from openmm import app
import mdtraj as md
import matplotlib.pyplot as plt

def add_residue_and_atoms(top,residue:str,chain):
    
    top.addResidue(residue,chain)
    residue = list(top.residues())[-1]
    atom1 = top.addAtom('CA',app.Element.getBySymbol('C'),residue)
    return atom1

def build_topology_ang(Chain_name,residue:list):
    """
       build a topology of angle
    """
    top = app.Topology()
    top.addChain(Chain_name)
    backbone = top._chains[-1]
    for i_resi in residue:
        add_residue_and_atoms(top,i_resi,backbone)
    atoms = list(top.atoms())
    for i in range(len(atoms)-1):
        top.addBond(atoms[i],atoms[i+1])
    return top

def setup_params(i_ang_idx,residue_params):
    """
       setup a pandas table format parameters for flexible local potential of bond angle
    """
    num_ang = len(i_ang_idx)
    pd_ang_idx_para_all = []
    if len(i_ang_idx[0]) == 3:
        columns = ['a1','a2','a3']
    elif len(i_ang_idx[0]) ==4:
        columns = ['a1','a2','a3','a4']
    for i in range(num_ang):
        pd_ang_idx = pd.DataFrame([i_ang_idx[i]],columns=columns)
        pd_ang_para = pd.DataFrame([residue_params])
        pd_ang_idx_para = pd.concat([pd_ang_idx,pd_ang_para],axis=1)
        pd_ang_idx_para_all.append(pd_ang_idx_para)
    pd_ang_idx_para_all = pd.concat(pd_ang_idx_para_all,axis=0)
    return pd_ang_idx_para_all

def construt_ang_simulation(midresidue,flp_bang_para,N_particles=3,T=300,box_len=5):
    """
       make a simulation configuration
    """
    system = mm.System()
    system.setDefaultPeriodicBoxVectors([box_len,0,0],[0,box_len,0],[0,0,box_len])
    for _ in range(N_particles):
        system.addParticle(137*unit.amu)
    # make topology
    top = build_topology_ang('ANG',['CGA',midresidue,'CGA'])
    # configure parameters
    atoms = list(top.atoms())
    ang_idxes = [[atoms[i].index,atoms[i+1].index,atoms[i+2].index] for i in range(len(atoms)-2)]
    pd_and_idx_para_all = setup_params(ang_idxes,flp_bang_para[midresidue])
    # set paramters 
    flex_loc_ang_force = flex_angle_term(pd_and_idx_para_all,force_group=0)
    system.addForce(flex_loc_ang_force)
    integrator = mm.LangevinIntegrator(T*unit.kelvin,1.0/unit.picosecond,2.0*unit.femtosecond)
    # create a simulatoin according topology, system and integrator
    simulation = app.Simulation(top,system,integrator)

    return simulation,system,top

def construt_dihd_simulation(midresidue_pair,flp_dihd_para,N_particles=4,T=300,box_len=5):
    """
       make a simulation configuration
    """
    system = mm.System()
    system.setDefaultPeriodicBoxVectors([box_len,0,0],[0,box_len,0],[0,0,box_len])
    for _ in range(N_particles):
        system.addParticle(137*unit.amu)
    # make topology
    top = build_topology_ang('DIHD',['CGA',midresidue_pair[0:3],midresidue_pair[3:],'CGA'])
    # configure parameters
    atoms = list(top.atoms())
    dihd_ang_idxes = [[atoms[i].index,atoms[i+1].index,atoms[i+2].index,atoms[i+3].index] for i in range(len(atoms)-3)]
    pd_dihdang_idx_para_all = setup_params(dihd_ang_idxes,flp_dihd_para)
    # set paramters 
    flex_loc_dihd_ang_force = flex_dihd_term(pd_dihdang_idx_para_all,force_group=0)
    system.addForce(flex_loc_dihd_ang_force)
    integrator = mm.LangevinIntegrator(T*unit.kelvin,1.0/unit.picosecond,2.0*unit.femtosecond)
    # create a simulatoin according topology, system and integrator
    simulation = app.Simulation(top,system,integrator)
    return simulation,system,top

# get energy and force for selected potential energy function
def get_energy_and_force(simulation,positions,groups_index):
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy=True,groups={groups_index})
    force = simulation.context.getState(getForces=True,groups={groups_index}).getForces(asNumpy=True)
    #force = [np.sqrt(f.dot(f)) for f in force]
    return state.getPotentialEnergy()._value,force

def compute_flex_loc(flex_ang_params_class,midresidue):

    bang_x = flex_ang_params_class.bond_ang_x
    bang_y = flex_ang_params_class.bond_ang_y
    bang_y2 = flex_ang_params_class.bond_ang_y2
    for i in range(len(bang_x)-1):
        theta_lo = bang_x[i]
        theta_hi = bang_x[i+1]
    
        y_lo = bang_y[midresidue][i]
        y_hi = bang_y[midresidue][i+1]
    
        y2_lo = bang_y2[midresidue][i]
        y2_hi = bang_y2[midresidue][i+1]
        if i == 0: 
           theta = np.arange(theta_lo,theta_hi,0.0005)
           y = flex_ang_params_class.cubic_spline(theta,theta_lo,theta_hi,y_lo,y_hi,y2_lo,y2_hi)
        else:
            i_theta = np.arange(theta_lo,theta_hi,0.0005)
            i_y = flex_ang_params_class.cubic_spline(i_theta,theta_lo,theta_hi,y_lo,y_hi,y2_lo,y2_hi)
            theta = np.hstack((theta,i_theta))
            y = np.hstack((y,i_y))
    return theta,y-np.min(y)

def generate_flex_loc_ang_fig(simulation,top,ang,midresidue,save_folder,flex_loc_potential_params):
    """
       generate the energy figure of a residue-dependent flexible local potential.
    """
    ang_all = []
    energy_all = []
    for iang in ang:
        x = 0 + np.cos(iang)
        z = 0 + np.sin(iang)
        positions = np.array([[0.38, 0, 0],[0,0,0],[x,0,z]])
        positions = positions * unit.nanometer
        energy,_ = get_energy_and_force(simulation,positions,0)
        traj = md.Trajectory(positions,top)
        ang_i = md.compute_angles(traj,[[0, 1, 2]])
        ang_all.append(ang_i[0][0])
        energy_all.append(energy)
    fig,ax = plt.subplots(figsize=(4.6,4))
    ax.plot(np.array(ang_all),np.array(energy_all)/_kcal_to_kj,'-',linewidth=1.5,color='k',alpha=1,label='OpenMM')
    cubicspine_theta,cubicspine_y = compute_flex_loc(flex_loc_potential_params,midresidue)
    ax.plot(cubicspine_theta,cubicspine_y,'-',linewidth=4,color='r',alpha=0.5,label='Python')
    ax.set_xlabel(r'$\theta$',fontsize=20)
    ax.set_ylabel('Energy (kcal/mol)',fontsize=20)
    plt.legend(fontsize=13)
    plt.savefig('%s/flp_ang/flp_ang_%s.png'%(save_folder,midresidue),dpi=600,bbox_inches='tight')
    plt.close()

def generate_flex_loc_dihd_fig(residue_pair,flex_ang_params,save_folder):
    flp_dihd_para = flex_ang_params.flp_bond_dihd_params
    dihd_ang = np.linspace(-4*np.pi,4*np.pi,400)
    energy_dihd = []
    dihd_all = []
    for i in range(len(dihd_ang)):
        ang = -dihd_ang[i]
        x = 0 + np.cos(ang)
        z = 0 + np.sin(ang)
        positions = np.array([[0.38,0,0],[0,0,0],[0,0.38,0],[x,0.38,z]])
        positions = positions * unit.nanometer
        simulation_dihd,system_dihd,top_dihd = construt_dihd_simulation(residue_pair,flp_dihd_para[residue_pair])
        ener, _ = get_energy_and_force(simulation_dihd,positions,0)
        traj = md.Trajectory(positions,top_dihd)
        dihd = md.compute_dihedrals(traj,[[0, 1, 2, 3]])
        energy_dihd.append(ener)
        dihd_all.append(dihd[0][0])
    y_dihd = flex_ang_params.flexi_dihd_energy(dihd_ang,flp_dihd_para[residue_pair]) 
    fig,ax = plt.subplots(figsize=(4.6,4))
    ax.plot(dihd_ang,np.array(energy_dihd)/_kcal_to_kj,'-',linewidth=1.5,color='k',alpha=1,label='OpenMM')
    ax.plot(dihd_ang,(y_dihd-np.min(y_dihd))/_kcal_to_kj,'-',linewidth=4,color='r',alpha=0.5,label='Python')
    ax.set_xlabel(r'$\theta$',fontsize=20)
    ax.set_ylabel('Energy (kcal/mol)',fontsize=20)
    plt.legend(fontsize=13)
    plt.savefig('%s/flp_dihd/flp_dihd_%s.png'%(save_folder,residue_pair),dpi=600,bbox_inches='tight')
    plt.close()
# Create a system for only angle
N_particles = 3
T = 300
box_len = 5 * unit.nanometer
_kcal_to_kj = 4.1840
save_folder = 'result'
# get parameter of flexible loc potential
flp_para_path = '../para/flexible_local_abeta.para'
flex_ang_params = parser_flp_para()
flex_ang_params.get_corr_flex_ang_para(flp_para_path)
flp_bang_para = flex_ang_params.flp_bond_ang_params
flp_dihd_para = flex_ang_params.flp_bond_dihd_params
residue_list = list(flp_bang_para.keys())

for resi_i in residue_list:
    simulation_ang,system_ang,top_ang = construt_ang_simulation(midresidue=resi_i,flp_bang_para=flp_bang_para)
    ang = np.linspace(1.1,3,100)
    generate_flex_loc_ang_fig(simulation_ang,top_ang,ang,resi_i,save_folder,flex_ang_params)

residue_pair = list(flp_dihd_para.keys())

for i_resi_pair in residue_pair:
    generate_flex_loc_dihd_fig(i_resi_pair,flex_ang_params,save_folder)
