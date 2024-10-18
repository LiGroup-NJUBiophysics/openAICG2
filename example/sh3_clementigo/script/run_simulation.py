import openmm as mm
from openmm import app
from openmm import unit
import sys
import numpy as np
import pandas as pd
import mdtraj as md
import json
import os 

from openaicg2.forcefield.aicgmodel import AICG2Model
from openaicg2 import utils 

T = 300
tot_simu_steps = 50000000
report_period = 1000
friction = 1
timestep = 20

pdb = md.load('../input/sh3_clementigo.pdb')
psf = md.load_psf('../input/sh3_clementigo.psf')

# complement the bond in topology
bonds = psf._bonds
top = pdb.topology.to_openmm()
rdtop = utils.RedefineTopology()
rdtop.redefine_bond(top,bonds)

# load native information
ParserNinfo=utils.ParserNinfo()
ParserNinfo.get_ninfo('../input/sh3_clementigo.ninfo')

model = AICG2Model()
model.create_system(top,use_pbc=True,
                    box_a=100, box_b=100, box_c=100,
                    nonbondedMethod=app.CutoffPeriodic,
                    remove_cmmotion=False)
# append native information to model
model.append_ff_params(ParserNinfo)
np.save('../output/native_contact.npy',model.protein_intra_contact[['a1','a2','sigma']].to_numpy())
# get exclusion for nonbonded interaction
model.get_exclusion(exclude_nat_con=True)
#print(len(model.extraexclusions))
# add force to system
model.add_protein_bond(force_group=0)
model.add_protein_harmonic_angle(force_group=1)
model.add_protein_native_dihedral(force_group=2)
model.add_protein_native_pair(force_group=3)
model.add_excluded(force_group=4)

# create a simulation
integrator = mm.LangevinIntegrator(T*unit.kelvin,friction/unit.picosecond,timestep*unit.femtosecond)
init_coord = pdb.xyz[0,:,:] * unit.nanometer
model.set_simulation(integrator, platform_name='CUDA',init_coord=init_coord)
model.move_COM_to_box_center(use_pbc=False)
model.simulation.context.setVelocitiesToTemperature(T*unit.kelvin)
model.simulation.minimizeEnergy()
model.add_reporters(tot_simu_steps, report_period, 
                    output_traj_name='../output/sh3_clementigo_%d'%T,report_traj_format='dcd'
                    ,report_traj=True,report_state_log=True)
print('Running simulation')
model.simulation.step(tot_simu_steps)
