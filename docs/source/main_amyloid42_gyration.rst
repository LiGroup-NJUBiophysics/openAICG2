Example 2. Simulating the Amyloid-:math:`\beta42` monomer.
===========================================================
Set simulation
--------------
This example will demonstrate how to create a two-bead coarse-grained model using the default AICG2+ force field of OpenAICG2. We choose the Amyloid $\beta 42$ monomer as the model protein. The Amyloid $\beta 42$ monomer is an IDP (intrinsically disordered protein), with most of its structure being unordered. However, when its monomer concentration exceeds the critical concentration, Amyloid $\beta 42$ will aggregate and deposit, forming relatively ordered structures. Next, we will use OpenAICG2 to run a Langevin dynamics simulation to calculate its radius of gyration.

.. code-block:: python

   # add required packages
   import openmm as mm
   from openmm import app
   from openmm import unit
   import mdtraj as md
   import json
   from openaicg2.forcefield.aicgmodel import AICG2Model
   from openaicg2 import utils

Load the parameters required for the simulation, including temperature, friction coefficient, time step, total simulation time, total number of steps, trajectory recording interval, output file name, etc.

.. code-block:: python

   # load simulation parameters from json file
   simu_params_path = '../input/simulationparams.json'
   configure_params = json.load(open(simu_params_path,'r'))
   T = configure_params['Temperature'] 
   friction = configure_params['friction']
   timestep = configure_params['timestep']
   tot_simu_steps = configure_params['total_mdsteps']
   report_period = configure_params['report_period']
   output_file = configure_params['output_file_name']
   platform_type = configure_params['platform_type']
   initial_pdb_path = configure_params['initial_pdb']
   monomer_psf_path = configure_params['monomer_psf']
   box = configure_params['box_vector']
   native_info_path = configure_params['native_information']
   lambdakh_scale = configure_params['lambdakh_scale']
   cutoff_kh = configure_params['cutoff_kh']

Initialize the topology from the PDB and PSF files (if the PDB file includes the CONNECT section, you do not need to load the PSF to complete the topology in the PDB).

.. code-block:: python 

   # load pdb and psf
   pdb = md.load(initial_pdb_path)
   psf = md.load_psf(monomer_psf_path)
   top = pdb.topology.to_openmm()
   top._bonds = []
   bonds = psf._bonds
   # Refine the bond in topology
   redefine_top = utils.RedefineTopology()
   redefine_top.redefine_bond(top,bonds)

Use the ParserNinfo in the utils to parse the native information file to obtain the force field parameters and store them in the ParserNinfo class variable.
# load parameter in native information file

.. code-block:: python 
   
   ParserNinfo=utils.ParserNinfo()
   ParserNinfo.get_ninfo(native_info_path)

After initialization, use the AICG2Model class to create an AICG model, then create an OpenMM system based on the model. The specific parameters include the initialized top, use_pbc=True to control whether nonbonded interactions are used, and nonbondedMethod to set the properties of nonbonded interactions, including whether they are truncated and whether periodic boundary conditions are used. After creating the system, we use the append_ff_params function to input the parsed force field parameters from ParserNinfo, and then call add_all_default_ener_function to add all the necessary default energy functions. The parameter oriented_Hbond set to False indicates that orientation-dependent hydrogen bonds are not added.

.. code-block:: python

   # create model
   model = AICG2Model()
   model.create_system(top,use_pbc=True,
                    box_a=box['x'], box_b=box['y'], box_c=box['z'],
                    nonbondedMethod=app.CutoffPeriodic,
                    remove_cmmotion=False)
   # input native information
   model.append_ff_params(ParserNinfo)
   # add all default energy function
   model.add_all_default_ener_function(oriented_Hbond=False,cutoffkh=cutoff_kh,
                                    kh_epsilon_scale=lambdakh_scale,temperature=T)

Create a Langevin integrator using OpenMM, and then build the simulation.

.. code-block:: python

   # create simulation
   integrator = mm.LangevinMiddleIntegrator(T*unit.kelvin,friction/unit.picosecond,timestep*unit.femtosecond)
   init_coord = pdb.xyz[0,:,:] * unit.nanometer
   model.set_simulation(integrator, platform_name=platform_type,properties={'Precision': 'mixed'},init_coord=init_coord)
   model.move_COM_to_box_center(use_pbc=False)
   model.simulation.context.setVelocitiesToTemperature(T*unit.kelvin)
   model.simulation.minimizeEnergy()
   # reporter trajectory and log about simulation
   model.add_reporters(tot_simu_steps, report_period,
                    output_traj_name='../output/%s'%(output_file),report_traj_format='dcd'
                    ,report_traj=True,report_state_log=True)

.. code-block:: python 

   print('Running simulation!!!')
   try:
      model.simulation.step(tot_simu_steps)
   except:
      print('simulation interruption')
      model.save_system('../output/system_%s.xml'%(output_file))
      model.save_state('../output/state_%s.xml'%(output_file))

Compute the radius of gyration :math:`R_{g}`
--------------------------------------------

The standard polymer physics formula for calculating the radius of gyration is as follows:

.. math::

   R_g = \sqrt{\frac{1}{N}<\sum_{i=1}^{N}(r_i - r_{CM})>}


In the equation, $N$ refers to the number of particles, :math: `r_{i}` is the coordinate of particle :math: `i`, :math: `r_{CM}`  is the center of mass coordinate, and `<...>` denotes ensemble average.

.. code-block:: python 

   import numpy as np
   import matplotlib.pyplot as plt

   _A_to_nm = 0.1
   lambdascale = 1.13
   traj = md.load('../output/monomer.dcd',top='../output/monomer.pdb')
   rg = md.compute_rg(traj,masses=None)
   Rg = np.mean(rg)/_A_to_nm
   print('radius of gyration:',Rg,'angstrom')
