import sys

from openaicg2.forcefield.functionterms.nonbonded_terms import oriented_dependent_Hbond_term

import openmm as mm
import numpy as np
import pandas as pd
from openmm import unit
from openmm import app

# get energy and force for selected potential energy function
def get_energy_and_force(simulation,positions,groups_index):
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy=True,groups={groups_index})
    force = simulation.context.getState(getForces=True,groups={groups_index}).getForces(asNumpy=True)
    #force = [np.sqrt(f.dot(f)) for f in force]
    return state.getPotentialEnergy()._value,force

def oriented_dependent_Hbond_term2(odHbond_info,use_pbc=False,cutoff=2.5,force_group=1):
    """a orientation dependent Lenard-Jones type potential for hydrogen bonds

    Parameter
    ---------
    odHbond_info: pd.DataFrame
        Information for all native hydrogen bonds are expressed as a N*Nine tabel.
        N row are represent the amount of Hbond. 9 columns are index of a1(atom1), 
        a2, a3, a4, a5, a6 and epsilon, sigma, cutoff(truncation distance/sigma), respectively

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default), then pbc would not apply to force.

    cutoff: float
        cutoff = truncation distance / rest length. If cutoff is None, the potential will not be trancated.

    force_group: int
        Force group.
    
    Return
    ------
    oriented_dependent_Hbond_force: Force
        OpenMM Force object
    """
    # energy function of Hbond
    energy_func = f"""  A;A = (1/(1+10*s1_quartic) + 1/(1+10*s2_quartic))/2;
                        s1_quartic= (1-c1)*(1+c1)*(1-c1)*(1+c1);
                        s2_quartic= (1-c2)*(1+c2)*(1-c2)*(1+c2);  
                        c1=(dotxv1 + dotyv1 + dotzv1)/(disp1*distance(p2,p5));
                        c2=(dotxv2 + dotyv2 + dotzv2)/(disp2*distance(p5,p2));
                        dotxv1=px1*dx25;dotyv1=py1*dy25;dotzv1=pz1*dz25;
                        dotxv2=px2*dx52;dotyv2=py2*dy52;dotzv2=pz2*dz52;
                        disp1= sqrt(px1^2 + py1^2 + pz1^2);
                        disp2= sqrt(px2^2 + py2^2 + pz2^2);
                        px1=(d1y*d2z-d1z*d2y);py1=(d1z*d2x-d1x*d2z);pz1=(d1x*d2y-d1y*d2x);
                        px2=(d3y*d4z-d3z*d4y);py2=(d3z*d4x-d3x*d4z);pz2=(d3x*d4y-d3y*d4x);
                        d1x=x1-x2;d1y=y1-y2;d1z=z1-z2;d2x=x3-x2;d2y=y3-y2;d2z=z3-z2;
                        d3x=x4-x5;d3y=y4-y5;d3z=z4-z5;d4x=x6-x5;d4y=y6-y5;d4z=z6-z5;
                        dx25=x5-x2;dy25=y5-y2;dz25=z5-z2;
                        dx52=x2-x5;dy52=y2-y5;dz52=z2-z5;"""
    energy_func = f"""step(distance(p5,p2)-sigma)*LJ*A+(1-delta(sigma-distance(p5,p2)))*step(sigma-distance(p5,p2))*(LJ+(1-A)*epsilon);
                                              LJ=epsilon*(5*(sigma/distance(p5,p2))^12-6*(sigma/distance(p5,p2))^10);
                                        A=(1/(1+10*sin(theta1)^4)+1/(1+10*sin(theta2)^4))/2;  
                                        theta1=pointangle(cx1,cy1,cz1,x2,y2,z2,x5,y5,z5);
                                        theta2=pointangle(cx2,cy2,cz2,x5,y5,z5,x2,y2,z2);
                                        cx1=px1+x2;cy1=py1+y2;cz1=pz1+z2;
                                        cx2=px2+x5;cy2=py2+y5;cz2=pz1+z5;
                                        px1=(d1y*d2z-d1z*d2y);py1=(d1z*d2x-d1x*d2z);pz1=(d1x*d2y-d1y*d2x);
                                        px2=(d3y*d4z-d3z*d4y);py2=(d3z*d4x-d3x*d4z);pz2=(d3x*d4y-d3y*d4x);
                                        d1x=x1-x2;d1y=y1-y2;d1z=z1-z2;d2x=x3-x2;d2y=y3-y2;d2z=z3-z2;
                                        d3x=x4-x5;d3y=y4-y5;d3z=z4-z5;d4x=x6-x5;d4y=y6-y5;d4z=z6-z5;"""
    if cutoff != None:
        energy_func = f"""step({cutoff}*sigma - distance(p5,p2))*""" + energy_func

    oriented_dependent_Hbond_force = mm.CustomCompoundBondForce(6,energy_func)
    oriented_dependent_Hbond_force.addPerBondParameter("epsilon")
    oriented_dependent_Hbond_force.addPerBondParameter("sigma")
    for row in odHbond_info.itertuples():
        oriented_dependent_Hbond_force.addBond([row.a2,row.a3 ,row.a1,row.a5,row.a6,row.a4], [row.epsilon, row.sigma])
    oriented_dependent_Hbond_force.setUsesPeriodicBoundaryConditions(use_pbc)
    oriented_dependent_Hbond_force.setForceGroup(force_group)
    return oriented_dependent_Hbond_force

# Create a simple system
box_len = 5*unit.nanometer
N_particle = 6
T = 300
system = mm.System()
system.setDefaultPeriodicBoxVectors([box_len,0,0],[0,box_len,0],[0,0,box_len])
for _ in range(N_particle):
    system.addParticle(137*unit.amu)
#positions = np.array([[2.88,2.5,2.5],[2.5,2.5,2.5],[2.5,2.88,2.5],[2.88,2.5,3.0],[2.5,2.5,3.0],[2.5,2.88,3.0]])
positions = np.array([[0.38,0,0],[0,0,0],[0,0.38,0],[0.38,0,0.9],[0,0,0.5],[0,0.38,0.9]])
positions = positions * unit.nanometer

# make topology 
num_residue = 2
num_atom_per_resi = 3
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

# initialize parameter of oriented dependent Hbond force and unit        
epsilon_hb = 1.58 # kcal/mol
sigma_hb = 5 # angstorm
cutoff = 3 # truncation distance/sigma
_A_to_nm_ = 0.1
_kcal_to_kj_ = 4.1840

# organize the hbond information for oriented dependent Hbond 
num_hbond = 1
pd_hb_info_all = []
for i in range(num_hbond):
    # Hbond involved indexes: a1 a2 a3 a4 a5 a6
    hbond_index = np.array([[0, 1, 2, 3, 4, 5,1,4]])
    # parameters of Hbond: epsilon, sigma, cutoff distance
    hbond_para = np.array([[epsilon_hb*_kcal_to_kj_,sigma_hb*_A_to_nm_,sigma_hb*_A_to_nm_*cutoff]])
    pd_hb_idx = pd.DataFrame(hbond_index,columns=['a1','a2','a3','a4','a5','a6','a7','a8'])
    pd_hb_para = pd.DataFrame(hbond_para,columns=['epsilon','sigma','cutoff']) 
    pd_hbond_info = pd.concat([pd_hb_idx,pd_hb_para],axis=1)
    pd_hb_info_all.append(pd_hbond_info)
# all information about Hbond
oriented_hbond_info = pd.concat(pd_hb_info_all,axis=0)

# create a orienteed dependent Hbond force 
Hbond_force = oriented_dependent_Hbond_term(oriented_hbond_info,use_pbc=True,force_group=0)
Hbond_force2 = oriented_dependent_Hbond_term2(oriented_hbond_info,use_pbc=True,force_group=1)
# add Hbond force to system
system.addForce(Hbond_force)
system.addForce(Hbond_force2)
integrator = mm.LangevinIntegrator(T*unit.kelvin,1.0/unit.picosecond,2.0*unit.femtosecond)
# create a simulatoin according topology, system and integrator
simulation = app.Simulation(top,system,integrator)

# 
position_parallel_not_plane = positions
ener_parallel_not_plane,force_parallel_not_plane = get_energy_and_force(simulation,position_parallel_not_plane,0)
ener_parallel_not_plane2,force_parallel_not_plane2 = get_energy_and_force(simulation,position_parallel_not_plane,1)
print(ener_parallel_not_plane,ener_parallel_not_plane2)
sys.exit()
position_parallel_in_plane = np.array([[0.38,0,0],[0,0,0],[0,0.38,0],[-0.88,0,0],[-0.5,0,0],[-0.5,0.38,0]])
ener_parallel_in_plane,force_parallel_in_plane = get_energy_and_force(simulation,position_parallel_in_plane,0)

# check the Hbond energy function in terms of energy 
try:
   assert(ener_parallel_not_plane==ener_parallel_in_plane*11)
except:
    print('Calculated hydrogen bond energies are wrong.')
