import pandas as pd
import openmm as mm
from openmm import unit
import numpy as np


def native_dihd_term(nat_dihd_info,use_pbc=False,force_group=1):
    """
    ordinary dihedral angle potential.

    Parameters
    ----------
    nat_dihd_info : pd.DataFrame
        Information for all native dihedral.     
    
    use_pbc : bool, optional 
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
    
    force_group : int
        Force group.
    """

    # create a custom Tortion force
    nat_tortion_force = mm.CustomTorsionForce("k_dih1*(1 - cos(dth)) + k_dih3*(1 - cos(3*dth));dth=theta-theta0")
    nat_tortion_force.addPerTorsionParameter("k_dih1")
    nat_tortion_force.addPerTorsionParameter("theta0")
    nat_tortion_force.addPerTorsionParameter("k_dih3")
    for row in nat_dihd_info.itertuples():
        nat_tortion_force.addTorsion(row.a1, row.a2, row.a3, row.a4, [row.k_dihd1,row.natdihd,row.k_dihd3])
    nat_tortion_force.setUsesPeriodicBoundaryConditions(use_pbc)
    nat_tortion_force.setForceGroup(force_group)

    return nat_tortion_force    

def aicg_dihd_term(aicg_dihd_info,use_pbc=False,force_group=1):
    """
    A structure-based local contact potential to describe the chirality of local interactions (dihedral angle of backbone) and reference from 
    "Li, W., Wang, W., & Takada, S. (2014). Energy landscape views for interplays among folding, binding, and allostery of calmodulin domains. 
    Proceedings of the National Academy of Sciences, 111(29), 10550-10555".

    Parameters
    ----------
    aicg_dihd_info : pd.DataFrame
        Information for all aicg dihedral angle
    
    use_pbc : bool, optional 
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
    
    force_group : int
        Force group.
    
    Return
    ------
    aicg_dihd_force : Force
        OpenMM Force object
    """
    energy_function = f'''-epsilon*exp(gex);
            gex=-dt_periodic^2/(2*width^2);
            dt_periodic = dt-floor((dt+{np.pi})/(2*{np.pi}))*(2*{np.pi});
            dt=theta - theta0
            '''    
    aicg_dihd_force = mm.CustomTorsionForce(energy_function)
    aicg_dihd_force.addPerTorsionParameter("epsilon")
    aicg_dihd_force.addPerTorsionParameter("theta0")
    aicg_dihd_force.addPerTorsionParameter("width")
    for row in aicg_dihd_info.itertuples():
        aicg_dihd_force.addTorsion(row.a1, row.a2, row.a3, row.a4, [row.epsilon, row.natdihd, row.width])
    aicg_dihd_force.setUsesPeriodicBoundaryConditions(use_pbc)
    aicg_dihd_force.setForceGroup(force_group)
    return aicg_dihd_force

def flex_dihd_term(flex_dihd_info,use_pbc=False,force_group=1):
    """a sequence dependence flexible local dihedral interaction 
    flex_dihd_info : pd.DataFrame
    Information for all aicg dihedral angle
    
    Parameters
    ----------
    flex_dihd_info : pd.DataFrame
        Information for all flexible local dihedral angle
        
    use_pbc : bool, optional 
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
    
    force_group : int
        Force group.
    
    Return
    ------
    flex_dihd_force : Force
        OpenMM Force object
    """
    flex_dihd_force = mm.CustomTorsionForce("C+c1*cos(theta)+s1*sin(theta)+ \
                                          c2*cos(2*theta)+s2*sin(2*theta)+\
                                          c3*cos(3*theta)+s3*sin(3*theta)-fdih_ener_corr")
    flex_dihd_force.addPerTorsionParameter("C")
    flex_dihd_force.addPerTorsionParameter("c1")
    flex_dihd_force.addPerTorsionParameter("s1")
    flex_dihd_force.addPerTorsionParameter("c2")
    flex_dihd_force.addPerTorsionParameter("s2")
    flex_dihd_force.addPerTorsionParameter("c3")
    flex_dihd_force.addPerTorsionParameter("s3")
    flex_dihd_force.addPerTorsionParameter("fdih_ener_corr")
    for _,row in flex_dihd_info.iterrows():
        a1 = row['a1']
        a2 = row['a2']
        a3 = row['a3']
        a4 = row['a4']
        para = row[4:]
        flex_dihd_force.addTorsion(a1,a2,a3,a4,para)
    flex_dihd_force.setUsesPeriodicBoundaryConditions(use_pbc)
    flex_dihd_force.setForceGroup(force_group)
    return flex_dihd_force
