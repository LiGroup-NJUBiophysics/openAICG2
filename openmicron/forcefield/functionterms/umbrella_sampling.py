import openmm as mm
from openmm import unit
import numpy as np
import pandas as pd

def umbrella_sampling_contact(cv_info,kums,q0,beta=50,lamda=1.2,force_group=12):
    """
    Implementing umbrella sampling for contacts

    Parameters
    ----------
    cv_info: pd.DataFrame
        The information for the contacts is expressed in a table format of N rows, 
        where N represents the number of contacts. Each row consists of three columns: 
        index (a1, a2), and the corresponding rest length (sigma) for each contact.
    
    kums: float
        The spring constant

    q0: float
        The center of each harmonic biasing potential. The values of q0 cover the range 
        from 0 to 1.

    beta: float
        The smoothing parameter that default value is 50 and unit is nm. 
    
    lamda: float
        The scale parameter of rest length. default value is 1.2.
    
    force_group: int
        Force group.
   
    Return
    ------
    ums_contact: openmm force
    """
    num_contact = cv_info.shape[0]
    cv_atom_idx_set = list(np.unique(cv_info[['a1','a2']].to_numpy().flatten()))
    for i,row in cv_info.iterrows():
        a1 = row.a1
        a2 = row.a2
        sigma = row.sigma
        a1_idx = cv_atom_idx_set.index(a1)+ 1
        a2_idx = cv_atom_idx_set.index(a2) + 1
        if i == 0:
            i_cv_dis = f'r{i}=distance(p{a1_idx},p{a2_idx});'
            i_cv_disrange = f"""rrange{i}=step(1+{lamda*sigma}-r{i});"""
            i_cv_ddis = f"""dr{i}=rrange{i}*(r{i}-{lamda*sigma});"""
            i_cv_qi = f"""contact{i}=rrange{i}/(1+exp({beta}*dr{i}));"""#step(1-r{a1_idx}{a2_idx}+{lambdar*h_cut})
            contact = f"""0.5*kums*(q-q0)^2;q={1/num_contact}*(contact{i}""" #{1/num_inter_contact_idx};0.5*kums*(q-q0)^2
        else:
            i_cv_dis += f'r{i}=distance(p{a1_idx},p{a2_idx});'
            i_cv_disrange += f"""rrange{i}=step(1+{lamda*sigma}-r{i});"""
            i_cv_ddis += f"""dr{i}=rrange{i}*(r{i}-{lamda*sigma});"""
            i_cv_qi += f"""contact{i}=rrange{i}/(1+exp({beta}*dr{i}));"""#step(1-r{a1_idx}{a2_idx}+{lambdar*h_cut})
            contact +=  f"""+contact{i}"""
    contact_function =  contact +');' +i_cv_qi + i_cv_ddis + i_cv_disrange + i_cv_dis
    ums_contact = mm.CustomCompoundBondForce(len(cv_atom_idx_set),contact_function)
    #ums_contact.addGlobalParameter('kums',kums)
    ums_contact.addPerBondParameter('kums')
    ums_contact.addPerBondParameter('q0')
    ums_contact.addBond(cv_atom_idx_set,[kums,q0])
    ums_contact.setForceGroup(force_group)
    return ums_contact
