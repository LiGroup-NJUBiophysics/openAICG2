import openmm as mm
import numpy as np
from openmm import app
from openmm import unit
import sys
import os
import pandas as pd

from openaicg2.utils import RedefineTopology
from openaicg2 import utils

__location__ = os.path.dirname(os.path.abspath(__file__))
_A_to_nm = 0.1
class SimulationSystem(object):
    '''
    Attributes
    -----------    

    '''
    def __init__(self):
        """
        Initialize
        """
        self.atoms = None
        
    def append_ff_params(self, ff_params,verbose=False):
        """
        The method can append new molecules by concatenating atoms and bonded interaction information saved in dataframes. 
        Reference from https://github.com/ZhangGroup-MITChemistry/OpenABC/blob/main/openabc/forcefields/cg_model.py

        Parameters
        ----------
        ff_params: force field parameters
            The object of a force field paramters including interaction information.
        
        Verbose: bool
            Whether to report the appended attributes.
        """
        chains = list(self.top._chains)
        mini_component_set = RedefineTopology()
        mini_component_set.get_mini_component_set(chains)
        num_mini_component_set = mini_component_set.num_mini_component_set
        num_atom_per_set = mini_component_set.num_atom_per_set
        # configure bond parameters 
        for i_set in range(num_mini_component_set):
            for each_attr_name in self.bonded_attr_name:
                if verbose:
                    print(f'Append attribute: {each_attr_name}')
                if hasattr(ff_params,each_attr_name):
                    if getattr(ff_params, each_attr_name) is not None:
                        new_attr = getattr(ff_params, each_attr_name).copy()
                        # reset the indices of atoms.
                        for each_col in ['a1','a2','a3','a4']:
                            if each_col in new_attr.columns:
                                new_attr[each_col] += i_set * num_atom_per_set
                        # set or conact index and parameters of molecular and interaction.
                        if hasattr(self, each_attr_name):
                            if getattr(self,each_attr_name) is None:
                                setattr(self, each_attr_name, new_attr)
                            else:
                                combined_attr = pd.concat([getattr(self, each_attr_name).copy(), new_attr],
                                                        ignore_index=True)
                                setattr(self, each_attr_name, combined_attr)
                        else:
                            setattr(self, each_attr_name,new_attr)
        # configure nonbonded parameters                    
        for i_set in range(num_mini_component_set-1):
            for j_set in range(i_set+1,num_mini_component_set):
                for each_attr_name in self.nonbonded_attr_name:
                    if verbose:
                        print(f'Append attribute: {each_attr_name}')
                    if hasattr(ff_params, each_attr_name):
                        #print(getattr(ff_params,each_attr_name),each_attr_name)
                        if getattr(ff_params,each_attr_name) is not None:
                            new_attr = getattr(ff_params, each_attr_name).copy()
                            new_attr.loc[:,['a1']] = new_attr.loc[:,['a1']] + i_set * num_atom_per_set  
                            new_attr.loc[:,['a2']] = new_attr.loc[:,['a2']] + j_set * num_atom_per_set
                        if hasattr(self, each_attr_name):
                            if getattr(self, each_attr_name) is None:
                                setattr(self, each_attr_name, new_attr)
                            else:
                                combined_attr = pd.concat([getattr(self, each_attr_name).copy(), new_attr],
                                                          ignore_index=True)
                                setattr(self, each_attr_name, combined_attr)
                        else:
                            setattr(self,each_attr_name,new_attr) 
                        
    def create_system(self,top,forcefield_template=None,use_pbc=True,box_a=100, box_b=100, box_c=100,nonbondedMethod=app.CutoffNonPeriodic,remove_cmmotion=False):
        """
        Create OpenMM system for simulation.
        Need to further add forces to this OpenMM system.

        Parameters
        ----------
        top: OpenMM Topology
            The OpenMM topology.

        forcefield_template: string 
            The path of OpenMM force fild xml file that define 
            the atom types and residue templates of system.
        
        use_pbc: bool
            Wheter to use periodic boundary condition (PBC).
        
        box_a: float 
            The length of the box along the x-axis is measured in angstroms.
        
        box_b: float
            The length of the box along the y-axis is measured in angstroms.
        
        box_c: float
            The length of the box along the z-axis is measured in angstroms.
        
        nonbondedMethod: Openmm method
            Set the method used for handling long range nonbonded interactions.
            Allowed values are NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.

        remove_commotion: bool
            Whether to remove center fo mass motions

        """
        forcefield_template = f'{__location__}/amyloid.xml'
        self.top = top
        self.atoms = list(self.top.atoms())
        self.chains = list(self.top.chains())
        self.use_pbc = use_pbc
        if self.use_pbc:
            box_vec = np.array([[box_a,0,0], [0,box_b,0], [0,0,box_c]])*_A_to_nm*unit.nanometer
            self.top.setPeriodicBoxVectors(box_vec)
        forcefield = app.ForceField(forcefield_template)
        self.system = forcefield.createSystem(self.top,nonbondedMethod=nonbondedMethod,
                                              nonbondedCutoff=2*unit.nanometer,
                                              removeCMMotion=remove_cmmotion)
        # check coarse-grained form
        atom_name = [i.name for i in self.atoms]
        atom_type = list(set(atom_name))
        if (len(atom_type) == 2 and 'CA' in atom_type and 'CB' in atom_type): 
            self.cgmodel_type = "two-bead"
        elif (len(atom_type) == 1 and 'CA' in atom_type):
            self.cgmodel_type = "one-bead"
        else:
            print("Error: The coarse-grained form not satisfy the requirment that two bead per residue or one per residue.")
        
    def get_exclusion(self, res_idx_dis=2,exclude_nat_con=True):
        """
        To get the exclusion that exclude nonbonded interactions when the distance of residue index is less than or equal to threshold value(res_idx_dis).

        Parameters
        ----------
        res_idx_dis: int
            The threshold distance of residue index and default value is 2。 
        
        exclude_nat_con: bool, optional
            If set to True (default), those atom pairs involved in native contacts, such as those formed by Go or hydrogen bonds,
            will be excluded from other LJ-type nonbonded potential calculations.

        """
        print('get_exclusion')
        self.atoms = list(self.top.atoms())
        exclusions_CB_CB = [[ai.index,aj.index] for i, ai in enumerate(self.atoms) for j, aj in enumerate(self.atoms[i+1:])
                          if abs(ai.residue.index-aj.residue.index) <=res_idx_dis
                          and aj.name=="CB" and ai.name =="CB" and
                          ai.residue.chain.index == aj.residue.chain.index]
        exclusions_CA_other =  [[ai.index,aj.index] for i, ai in enumerate(self.atoms) for j, aj in enumerate(self.atoms[i+1:])
                          if abs(ai.residue.index-aj.residue.index) <=res_idx_dis 
                          and'CA' in [ai.name,ai.name] and
                          ai.residue.chain.index == aj.residue.chain.index]
        self.exclusions = exclusions_CB_CB + exclusions_CA_other 
        if exclude_nat_con:
            self.extraexclusions = []
            if hasattr(self,'protein_intra_contact'):
                protein_intra_con_idx = self.protein_intra_contact.loc[:,['a1','a2']]
                self.extraexclusions += protein_intra_con_idx.values.tolist()
            if hasattr(self,'protein_inter_contact') :
                protein_inter_con_idx = self.protein_inter_contact.loc[:,['a1','a2']]
                self.extraexclusions += protein_inter_con_idx.values.tolist() 
            if hasattr(self, 'oriented_Hbond'):
                hbond_index = self.oriented_Hbond.loc[:,['a1','a2']]
                self.extraexclusions += hbond_index.values.tolist()
           
    def set_simulation(self, integrator, platform_name='CPU',properties={'Precision': 'mixed'},init_coord=None):
        """
        Set OpenMM simulation

        Parameters
        ----------
        integrator: Openmm Integrator
            OpenMM integrator.

        platform_name: str
            OpenMM simulation platform name. The available platforms are Reference or CPU or CUDA or OpenCL.
        
        properties: dict
            OpenMM simulation platform properties.
        
        init_coord: None or array-like
            Initial coordinate of system.
        """
        platform = mm.Platform.getPlatformByName(platform_name)
        print(f'Use platform: {platform_name}')
        if platform_name in ['CUDA', 'OpenCL']:
            if 'Precision' not in properties:
                properties['Precision'] = 'mixed'
                properties['DeviceIndex'] = '0'
            precision = properties['Precision']
            print(f'Use precision: {precision}')
            self.simulation = app.Simulation(self.top,self.system,integrator,platform,properties)
        else:
            self.simulation = app.Simulation(self.top, self.system, integrator, platform)
        if init_coord is not None:
            self.simulation.context.setPositions(init_coord)

    def move_COM_to_box_center(self,use_pbc=None):
        """
        Move center of mass (COM) to box center.
        """  
        print('Move center of mass (COM) to box center.')
        if use_pbc == None:
           state = self.simulation.context.getState(getPositions=True,enforcePeriodicBox=self.use_pbc)
        else:
           state = self.simulation.context.getState(getPositions=True,enforcePeriodicBox=use_pbc)
        positions = np.array(state.getPositions().value_in_unit(unit.nanometer))
        n_atoms = positions.shape[0]
        mass = [self.system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(n_atoms)]
        mass = np.array(mass)
        weights = mass/np.sum(mass)
        box_vec_a, box_vec_b, box_vec_c = self.system.getDefaultPeriodicBoxVectors()
        box_vec_a = np.array(box_vec_a.value_in_unit(unit.nanometer))
        box_vec_b = np.array(box_vec_b.value_in_unit(unit.nanometer))
        box_vec_c = np.array(box_vec_c.value_in_unit(unit.nanometer))
        box_center = 0.5*(box_vec_a + box_vec_b + box_vec_c)

        center_of_mass = np.average(positions, axis=0,weights=weights)
        positions = positions - center_of_mass + box_center
        self.simulation.context.setPositions(positions*unit.nanometer)

    def save_system(self, system_xml="system.xml"):
        """
        Save sytem in a readabel XML format.

        Parameters
        ----------
        system_xml: str
            Output path for system xml file 

        """
        with open(system_xml,'w') as output_writer:
            output_writer.write(mm.XmlSerializer.serialize(self.system))
    
    def save_state(self, state_xml="state.xml"):
        """
        Save state in a reabable XML format

        Parameters
        ----------
        state_xml: str
            Output path for state xml file
        
        """
        with open (state_xml,'w') as output_writer:
            state = self.simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, 
                                                     getEnergy=True, getParameters=True,
                                                     enforcePeriodicBox=self.use_pbc)
            output_writer.write(mm.XmlSerializer.serialize(state))
    
    def add_reporters(self, tot_simu_steps, report_period, output_traj_name='output',report_traj_format='dcd',report_traj=True,report_state_log=True):
        """
        Add reporters that produce trajectory in dcd or xtc format and state log for OpenMM simulation.

        Parameters
        ----------
        tot_sim_steps: int
            Total time steps of simulation.

        report_period: int
            Report period for trajectory and state log files.

        output_traj_name: str
            Output path and name for trajectory and state log files.
        
        report_traj_format: str
            Trajectory file format, which one can choose dcd ro xtc format.
        
        report_traj: bool
            Whether to output trajectory file.
        
        report_state_log: bool
            Whether to output sate log file.
            
        """
        positions = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        app.PDBFile.writeFile(self.top,positions,open('%s.pdb'%(output_traj_name),'w'))
        utils.gene_psf(output_traj_name,self.top)
        if report_traj_format == 'dcd' and report_traj:
           self.simulation.reporters.append(app.DCDReporter('%s.dcd'%output_traj_name,report_period,append=False))
        elif report_traj_format == 'xtc' and report_traj:
           self.simulation.reporters.append(app.XTCReporter('%s.xtc'%output_traj_name,report_period,append=False))
        if report_state_log:
           self.simulation.reporters.append(app.StateDataReporter('%s.log'%output_traj_name,report_period,step=True,
                            totalEnergy=True,potentialEnergy=True,kineticEnergy=True, 
                            temperature=True, progress=True,remainingTime=True,    
                            speed=True,totalSteps=tot_simu_steps,separator='\t',append=False))
           
    def auto_get_charged_atom(self):
        """
        automatically get the index and charge of charged atom.  
        """
        info_residue_charge = {'ASP':-1,'GLU':-1,'LYS':1,'ARG':1,'HIS':0.5,'HIE':0.5}
        self.charged_atoms = []
        for i_atoms in self.atoms:
            resi_name = i_atoms.residue.name
            if i_atoms.name == 'CB' and resi_name in info_residue_charge.keys():
                self.charged_atoms.append([i_atoms.index,info_residue_charge[resi_name]])

    
    
