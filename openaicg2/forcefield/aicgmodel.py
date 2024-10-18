import numpy as np
import pandas as pd
from openmm import app
from openmm import unit
import sys
import os

from openaicg2.forcefield.simulationsystem import SimulationSystem
from openaicg2.forcefield import functionterms
from openaicg2 import utils

__location__ = os.path.dirname(os.path.abspath(__file__))
_A_to_nm = 0.1
_kcal_to_kj = 4.1840

class AICG2Model(SimulationSystem):
    """
    A class for AICG2+ model

    **Attributes:**

    bonded_attr_name:  list
        The categories of interactions within chains. Such as ['protein_bonds','protein_harmonic_angles','protein_aicg13_angles',
        'protein_native_dihd','protein_aicg_dihd','protein_intra_contact']

    nonbonded_attr_name: list
        The categories of interactions between chains.
        ['protein_inter_contact']
    
    flp_bond_ang_params: pd.DataFrame
        The parameters of flexible local potential regarding bond angle.
    
    flp_bond_dihd_params: pd.DataFrame
        The parameters of flexible local potential regarding dihedral. 
    
    protein_bonds: pd.DataFrame
        The parameters of harmonica bond interaction
    
    protein_harmonic_angles: pd.DataFrame
        The parameters of harmonica angle interaction
    
    protein_aicg13_angles: pd.DataFrame
        The parameters of aicg13 angles
    
    protein_native_dihd: pd.DataFrame
        The parameters of native dihedral
    
    protein_aicg_dihd: pd.DataFrame
        The parameters of aicg dihedral
    
    protein_intra_contact: pd.DataFrame
        The parameters of intra-contact
    
    protein_inter_contact: pd.DataFrame
        The parameters of inter-contact
    """
    def __init__(self):
        """
        Initialize.
        """
        self.bonded_attr_name = ['protein_bonds','protein_harmonic_angles','protein_aicg13_angles',
                                'protein_native_dihd','protein_aicg_dihd','protein_intra_contact']
        self.nonbonded_attr_name = ['protein_inter_contact']
        
        
        self.parser_flp = utils.parser_flp_para()
        self.parser_flp.get_corr_flex_ang_para(f'{__location__}/para/flexible_local_abeta.para')
        self.flp_bond_ang_params = self.parser_flp.flp_bond_ang_params
        self.flp_bond_dihd_params = self.parser_flp.flp_bond_dihd_params


    def add_protein_bond(self, force_group=1):
        """
        Add protein bonds.

        Parameters
        ----------

        force_group: int
            Force group.
        """

        if hasattr(self, 'protein_bonds'):
            print('Add protein bonds.')
            force = functionterms.harmonic_bond_term(self.protein_bonds, False,force_group)
            self.system.addForce(force)
    
    def add_protein_harmonic_angle(self,k_angle_scale=4,force_group=2):
        """
        Add protein harmonic angles.
        
        Parameters
        ----------

        force_group: int
            Force group.
        """

        if hasattr(self, 'protein_harmonic_angles'):
            print('Add protein angles.')
            force = functionterms.harmonic_angle_term( self.protein_harmonic_angles,False,force_group)
            self.system.addForce(force)

    def add_protein_aicg13_angle(self,force_group=3):
        """
        Add protein aicg13 angles involved CA atom that describe the local structure in backbone chain 

        Parameters
        ----------

        force_group: int
            Force group.
        """

        if hasattr(self, 'protein_aicg13_angles'):
            print('Add protien aicg13 angles.')
            aicg13_ang_params = self.protein_aicg13_angles
            ang_idx = []
            ang_params = []
            for row in aicg13_ang_params.itertuples():
                if row.epsilon!=0:
                    ang_idx.append([row.a1,row.a2,row.a3])
                    ang_params.append([row.epsilon,row.r0,row.width])
            ang_idx = np.array(ang_idx).astype(np.int64)
            ang_params = np.array(ang_params).astype(np.float64)
            pd_ang_idx = pd.DataFrame(ang_idx,columns=['a1','a2','a3'])
            pd_ang_params = pd.DataFrame(ang_params,columns=['epsilon','r0','width'])
            aicg13_ang_params_invo_ca = pd.concat([pd_ang_idx,pd_ang_params],axis=1)
            force = functionterms.aicg13_angle_term(aicg13_ang_params_invo_ca,force_group=force_group)
            self.system.addForce(force)
    
    def add_protein_native_dihedral(self,force_group=4):
        """
        add native protein dihedrals for side bead.

        Parameters
        ----------
        force_group: int
            Force group.
        """
        if hasattr(self, 'protein_native_dihd'):
            print('Add protein native dihedral angle.')
            nat_dihd_params = self.protein_native_dihd
            dihd_idx = []
            dihd_params = []
            for row in nat_dihd_params.itertuples():
                if row.k_dihd1 !=0 and row.k_dihd3!=0:
                    dihd_idx.append([row.a1,row.a2,row.a3,row.a4])
                    dihd_params.append([row.k_dihd1,row.natdihd,row.k_dihd3])
            dihd_idx = np.array(dihd_idx).astype(np.int64)
            dihd_params = np.array(dihd_params).astype(np.float64)
            pd_dihd_idx = pd.DataFrame(dihd_idx,columns=['a1','a2','a3','a4'])
            pd_dihd_params = pd.DataFrame(dihd_params,columns=['k_dihd1','natdihd','k_dihd3'])
            nat_dihd_params_invo_cb = pd.concat([pd_dihd_idx,pd_dihd_params],axis=1)            
            force = functionterms.native_dihd_term(nat_dihd_params_invo_cb,force_group=force_group)
            self.system.addForce(force)
    
    def add_protein_aicg_dihedral(self,force_group=5):
        """
        add protein aicg dihedral angle for backbone chain.

        Parameters
        ----------
        force_group: int
            Force group.
        """
        if hasattr(self,'protein_aicg_dihd'):
            print('Add protein aicg dihedral.')
            aicg14_dihd_params = self.protein_aicg_dihd
            dihd_idx = []
            dihd_params = []
            for row in aicg14_dihd_params.itertuples():
                if row.epsilon!=0: 
                    dihd_idx.append([row.a1,row.a2,row.a3,row.a4])
                    dihd_params.append([row.epsilon, row.natdihd,row.width])
            dihd_idx = np.array(dihd_idx).astype(np.int64)
            dihd_params = np.array(dihd_params).astype(np.float64)
            pd_dihd_idx = pd.DataFrame(dihd_idx,columns=['a1','a2','a3','a4'])
            pd_dihd_params = pd.DataFrame(dihd_params,columns=['epsilon','natdihd','width'])
            nat_dihd_params_invo_ca = pd.concat([pd_dihd_idx,pd_dihd_params],axis=1)   
            force = functionterms.aicg_dihd_term(nat_dihd_params_invo_ca, force_group=force_group)
            self.system.addForce(force)
    
    def add_flexible_loc_angle(self,force_group=6):
        """
        Add flexible local potential for angle are composed of CA atoms in backbone.

        Parameters
        ----------
        force_group: int
            Force group.
        """

        print('Add flexible local potential of angle of backbone.')
        flp_ang_idx = []
        flp_ang_params = []
        all_ca_atoms = [ai for ai in self.atoms if ai.name=='CA']
        
        for i in range(len(all_ca_atoms)-2):
            a1 = all_ca_atoms[i]
            a2 = all_ca_atoms[i+1]
            a3 = all_ca_atoms[i+2]
            chain_idx = set([a1.residue.chain.index,a2.residue.chain.index,a3.residue.chain.index])
            if len(list(chain_idx)) == 1:
                flp_ang_idx.append([a1.index,a2.index,a3.index])
                flp_ang_params.append(self.flp_bond_ang_params[a2.residue.name])
        flp_ang_idx = np.array(flp_ang_idx).astype(np.int64)
        flp_ang_params = np.array(flp_ang_params).astype(np.float64)
        pd_flp_ang_idx = pd.DataFrame(flp_ang_idx,columns=['a1','a2','a3'])
        pd_flp_ang_params = pd.DataFrame(flp_ang_params)
        flp_ang_idx_params = pd.concat([pd_flp_ang_idx,pd_flp_ang_params],axis=1)
        force = functionterms.flex_angle_term(flp_ang_idx_params,force_group=force_group)
        self.system.addForce(force)
            
    def add_flexible_loc_dihedral(self,force_group=7):
        """
        Add flexible local potential for dihedral angle are composed CA atoms in backbone.

        Parameters
        ----------
        force_group: int
            Force group.
        """

        print('Add flexible local potential of dihedral of backbone')
        flp_dihd_idx = []
        flp_dihd_params = []
        all_ca_atoms = [ai for ai in self.atoms if ai.name=='CA']
        for i in range(len(all_ca_atoms)-3):
            a1 = all_ca_atoms[i]
            a2 = all_ca_atoms[i+1]
            a3 = all_ca_atoms[i+2]
            a4 = all_ca_atoms[i+3]
            chain_idx = set([a1.residue.chain.index,a2.residue.chain.index,a3.residue.chain.index,a4.residue.chain.index])
            if len(list(chain_idx)) == 1:
                flp_dihd_idx.append([a1.index,a2.index,a3.index,a4.index])
                centre_residue_pair_name = a2.residue.name + a3.residue.name
                flp_dihd_params.append(self.flp_bond_dihd_params[centre_residue_pair_name])
        flp_dihd_idx = np.array(flp_dihd_idx).astype(np.int64)
        flp_dihd_params = np.array(flp_dihd_params).astype(np.float64)
        pd_flp_dihd_idx = pd.DataFrame(flp_dihd_idx,columns=['a1','a2','a3','a4'])
        pd_flp_dihd_params = pd.DataFrame(flp_dihd_params,columns=['C','c1','s1','c2','s2','c3','s3','fdih_ener_corr'])
        flp_dihd_idx_params = pd.concat([pd_flp_dihd_idx,pd_flp_dihd_params],axis=1)
        force = functionterms.flex_dihd_term(flp_dihd_idx_params,force_group=force_group)
        self.system.addForce(force)


    def add_protein_native_pair(self,cutoff=2.5,force_group=8):
        """
        Add native contact pair. 

        Parameters
        ----------
        force_group: int
            Force group.
        """
        if hasattr(self, 'protein_inter_contact') or hasattr(self,'protein_intra_contact'):
            if hasattr(self, 'protein_inter_contact') and hasattr(self,'protein_intra_contact'):
                print('Add protein inter and intra native contact')
                inter_contact = self.protein_inter_contact.copy()
                intra_contact = self.protein_intra_contact.copy()
                protein_nat_con = pd.concat([intra_contact,inter_contact],axis=0) 
            elif hasattr(self,'protein_inter_contact') and not hasattr(self,'protein_intra_contact'):
                print('Add protien inter native contact')
                protein_nat_con = self.protein_inter_contact.copy()
            elif not hasattr(self, 'protein_inter_contact') and hasattr(self,'protein_intra_contact'):
                print('Add protien intra native contact')
                protein_nat_con = self.protein_intra_contact.copy()
            force = functionterms.go_contact_term(protein_nat_con,self.use_pbc,cutoff,force_group)
            self.system.addForce(force)



    def add_kim_hummer(self,path=f'{__location__}/para/kh.para',kh_model_symbol='D',T=300,kh_epsilon_scale=1.3,cutoff=2.5, rad_scale = 0.85,force_group=10):
        """
        Add kim hummer potential between CB and CB which there is no native contact

        Parameters
        ----------
        path: str
            The path of kim-hummer parameters file.

        kh_model_symbol: str    
            The parameters represent symbols of the Kim-Hummer parameters model. The model consists of
            six types: A, B, C, D, E, and F. Specific numerical values are referenced from "Journal 
            of Molecular Biology, 2008, 375(5): 1416-1433."

        T: float
            Temperature

        kh_epsilon_scale: float
            The parameter scale the size of kim-hummer epsilon.

        cutoff: float
            cutoff = truncated distance / rest length

        rad_scale: float
            scaling factor to radii(sigma) and default value is 0.85
        
        force_group: int
            Force group.
        """
        print('Add kim-hummer force')

        resi_type,epsilon_KH,sigma_KH = utils.parser_kh_params(path,kh_model_symbol,kh_epsilon_scale,T)
        num_type_resi = len(resi_type.keys())
        sigma_KH = sigma_KH * rad_scale
        atom_resi_types = []
        for ai in self.atoms:
            atom_resi_types.append(resi_type[ai.residue.name])

        if hasattr(self,'extraexclusions'):
            extraexclusions = self.extraexclusions
        else:
            extraexclusions = None  
        force = functionterms.kim_hummer_term(atom_resi_types,epsilon_KH,sigma_KH,
                                              self.exclusions,extraexclusions,use_pbc=self.use_pbc,
                                              cutoff=cutoff,force_group=force_group) 
        self.system.addForce(force)
            
    def add_excluded(self, epsilon=0.2,sigma = 3.8,cutoff=2.5,rad_scale=1,force_group=11):
        """
        Add excluded interaction between atoms.

        Parameters
        ----------
        cutoff: float
            cutoff = truncated distance / rest length

        rad_scale: float
            scaling factor to radii(sigma) and default value is 0.85
        force_group: int
            Force group.
        """
        print('Add excluded force')
        atom_type = []
        exv_params = utils.parser_exv_params(f'{__location__}/para/exv.para')
        numbeadtype = len(exv_params.keys()) - 2
        keys = list(exv_params.keys())
        epsilon_map = np.zeros((numbeadtype,numbeadtype),dtype=float)
        sigma_map = np.zeros((numbeadtype,numbeadtype),dtype=float)
        epsilon_map[:] = exv_params['exv_coef']
        for i in range(numbeadtype):
            for j in range(numbeadtype):
                sigma_map[i,j] = (exv_params[keys[i]] + exv_params[keys[j]])/2
        for ai in self.atoms:
            atom_type.append(keys.index(ai.residue.name))   

        if hasattr(self,'extraexclusions'):
            extraexclusions = self.extraexclusions
        else:
            extraexclusions = None  
        sigma_map = sigma_map * rad_scale
        force = functionterms.excluded_term(atom_type,epsilon_map,sigma_map,
                                            self.exclusions,extraexclusions,use_pbc=self.use_pbc,
                                            cutoff=cutoff,force_group=force_group)
        self.system.addForce(force)
   

    def add_debye_huckel(self,dieletric_constant=80,ion_strength=0.02,temperature=300,extra_charged_atom=None,cutoff=20,force_group=13):
        """
        Add electrical potential.
        
        Parameters
        ----------
        ion_strength: float
            Salt concentration.
        
        temperature: float
            Temperature.

        extra_charged_atom: list
            This represents the additional atoms requiring charge to be added. It consists of N*2 lists, where the first column represents the atom index, 
            and the second column represents the atom's charge.
        
        cutoff: float
            cutoff = truncated distance / rest length
        
        force_group: int
            Force group. 

        """
        print('Add debye huckel potential')
        self.auto_get_charged_atom()
        if extra_charged_atom is not None:
            self.charged_atoms = self.charged_atoms + extra_charged_atom
        pair_index_charged_atoms = [[i_ind_q[0],j_ind_q[0],i_ind_q[1],j_ind_q[1]] for i,i_ind_q in enumerate(self.charged_atoms[:-1]) 
                                                                                for j,j_ind_q in enumerate(self.charged_atoms[i+1:])]
        force = functionterms.debye_Huckel_bond_form(pair_index_charged_atoms,self.exclusions,dieletric_constant=dieletric_constant,
                                                     ion_strength=ion_strength,T=temperature,
                                                     use_pbc=self.use_pbc,cutoff=cutoff,force_group=force_group)
        self.system.addForce(force)
    
    def add_all_default_ener_function(self,oriented_Hbond=False,cutoff_hbond=2.5,cutoff_go=2.5,cutoff_kh=2.5,cutoff_exv=2.0,kh_epsilon_scale=1.3,temperature=300):
        """
        Add all default energy function to create a aicg force field.

        Parameters
        ----------
        oriented_Hbond: bool
            Whether to add orientation-dependent hydrogen bond between CA atoms for force field. 
            If False, the hydrogen would not apply to f  rce field.  
        
        cutoff_hbond: float
            cutoff_hbond = truncated distance / rest length
        
        cutoff_go: float
            cutoff_go =  truncated distance / rest length
        
        cutoff_kh: float
            cutoff_kh =  truncated distance / rest length

        cutoff_exv: float
            cutoff_exv=  truncated distance / rest length

        kh_epsilon_scale: float
            The parameter scale the size of kim-hummer epsilon.

        """
        self.add_protein_bond(force_group=1)
        self.add_protein_harmonic_angle(force_group=2)
        self.add_protein_aicg13_angle(force_group=3)
        self.add_flexible_loc_angle(force_group=4)
        self.add_protein_native_dihedral(force_group=5)
        self.add_protein_aicg_dihedral(force_group=6)
        self.add_flexible_loc_dihedral(force_group=7)
        self.add_protein_native_pair(force_group=9) # cutoff=cutoffgo
        self.get_exclusion(exclude_nat_con=True)
        self.add_kim_hummer(cutoff=cutoff_kh,force_group=10)
        self.add_excluded(cutoff=cutoff_exv,force_group=11)
        atoms_per_chain = list(self.chains[0].atoms())
        terminal_charge = []
        for i in range(len(self.chains)):
            index_N = len(atoms_per_chain) * i + 0
            index_C = len(atoms_per_chain) * i + 76
            terminal_charge.append([index_N,1.0])
            terminal_charge.append([index_C,-1.0])
        self.add_debye_huckel(ion_strength=0.02, temperature=temperature,extra_charged_atom=terminal_charge,force_group=11)
    
