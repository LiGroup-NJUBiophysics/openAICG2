Example 1. protein sh3 based on one bead model
==============================================
To use this example on your local computer, you need to clone or download the folder example from the GitHub repository of openmicron.

.. code-block:: python
   
   # import required packages
   import openmm as mm
   from openmm import app
   from openmm import unit
   import numpy as np
   import pandas as pd
   import mdtraj as md
   from openmicron.forcefield.aicgmodel import AICG2Model
   from openmicron import utils

Firstly, initialize the topology from the PDB and PSF files (if the PDB file includes the CONNECT section, you do not need to load the PSF to complete the topology in the PDB). Use the ParserNinfo from utils to parse the native information file to obtain the force field parameters, and store them in the ParserNinfo class variable.

.. code-block:: python

   # Initialize
   T = 300
   tot_simu_steps = 5
   report_period = 1000
   friction = 1
   timestep = 20

   # load pdb and psf files
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

After initializing the topology and native information file, create an AICG model using the AICG2Model class, and then create an OpenMM system based on the model. The specific parameters include the initialized top, use_pbc to control whether nonbonded interactions are used, and nonbondedMethod to set the properties of nonbonded interactions, including whether they are truncated and whether periodic boundary conditions are used.

.. code-block:: python

   model = AICG2Model()
   model.create_system(top,use_pbc=True,
                    box_a=100, box_b=100, box_c=100,
                    nonbondedMethod=app.CutoffPeriodic,
                    remove_cmmotion=False)

After creating the system, we use the append_ff_params function to input the parsed force field parameters from ParserNinfo. For non-bonded interaction calculations, besides excluding atom pairs with residue indices 
, native contacts are also excluded, controlled by the exclude_nat_con parameter of the get_exclusion function. Then, add the various forces in sequence. It is important to assign a different force group to each force for easier subsequent force field analysis.

.. code-block:: python

   # append native information to model
   model.append_ff_params(ParserNinfo)
   np.save('../output/native_contact.npy',model.protein_intra_contact[['a1','a2','sigma']].to_numpy())
   # get exclusion for nonbonded interaction
   model.get_exclusion(exclude_nat_con=True)
   # add force to system
   model.add_protein_bond(force_group=0)
   model.add_protein_harmonic_angle(force_group=1)
   model.add_protein_native_dihedral(force_group=2)
   model.add_protein_native_pair(force_group=3)
   model.add_excluded(force_group=4)

Use OpenMM to create a Langevin integrator, then select the computing platform, such as CPU or CUDA, and build the OpenMM simulation.

.. code-block:: python

   # create a integrator
   integrator = mm.LangevinIntegrator(T*unit.kelvin,friction/unit.picosecond,timestep*unit.femtosecond)
   init_coord = pdb.xyz[0,:,:] * unit.nanometer
   # create a simulation
   model.set_simulation(integrator, platform_name='CPU',init_coord=init_coord)
   model.move_COM_to_box_center(use_pbc=False)
   model.simulation.context.setVelocitiesToTemperature(T*unit.kelvin)
   model.simulation.minimizeEnergy()

Call the `add_reporters` function inherited to `simulationsystem` to set the input and output. `tot_simu_steps` specifies the total number of simulation steps, `report_period` controls the interval at which information is recorded, `output_traj_name` sets the output file path and name, `report_traj_format` sets the output trajectory format, with options of "dcd" and "xtc". `report_traj` and `report_state_log` respectively set whether to output the trajectory file and the log file. Finally, run the simulation.

.. code-block:: python 

   model.add_reporters(tot_simu_steps, report_period,
                    output_traj_name='../output/sh3_clementigo',report_traj_format='dcd'
                    ,report_traj=True,report_state_log=True)
   print('Running simulation')
   model.simulation.step(tot_simu_steps)
