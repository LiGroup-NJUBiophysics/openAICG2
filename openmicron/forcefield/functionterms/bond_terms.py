import pandas as pd
import openmm as mm
from openmm import unit

def harmonic_bond_term(bond_info,use_pbc=False,force_group=1):
    """
    Harmonic bond term.

    Parameters
    ----------
    bond_info : pd.DataFrame
        Information for all bonds.
    
    use_pbc : bool, optional
        Whether to use periodic boundary conditions. If False (default),then pbc would not apply to Harmonic bond force. 

    force_group : int
        Force group 
    
    Return
    ------
    bond_force : Force 
        OpenMM Harmonic bond force object
    """

    harmonic_bond_force = mm.HarmonicBondForce()
    for row in bond_info.itertuples():
        atom1 = row.a1
        atom2 = row.a2
        r0 = row.r0
        k_bond = row.k
        harmonic_bond_force.addBond(atom1,atom2,r0,k_bond)
    harmonic_bond_force.setUsesPeriodicBoundaryConditions(use_pbc)
    harmonic_bond_force.setForceGroup(force_group)

    return harmonic_bond_force

