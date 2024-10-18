import openmm as mm
from openmm import unit
import numpy as np
import pandas as pd

_kcal_to_kj = 4.1840
_A_to_nm = 0.1
def fix_boundary(atoms,xbox,ybox,zbox,boxsigma=7.5,kbox=10,force_group=1):
    """
    Implementing elastic boundary conditions, which means making particles bounce back when they hit the wall.
    
    Parameters
    ----------
    atoms: list
        all atoms in system.
    xbox: float
        The box length in x-axis.

    ybox: float
        The box length in y-axis.

    zbox: float
        The box length in z-axis.

    boxsigma: float
        Sigma of box represents the sharpness of the wall. Default value is 7.5 \AA

    kbox: float
        The strength of wall. Defalut value is 10 kcal/mol

    force_group: int
        Force group.

    Return
    ------
    fix_boundary: openmm force
        A rectangular box of fixed dimensions.
    """
    
    energy_function = f"""ex+ey+ez;
                          ez=ez1 + ez2 - kbox*(1/0.8)^12-emin;
                          ez1=step(dz-0.8*boxsigma)*step(3*boxsigma-dz)*kbox*(boxsigma/dz)^12;
                          ez2=step(0.8*boxsigma-dz)*kbox*(boxsigma/(0.8*boxsigma))^12*(1 + 12*(0.8*boxsigma - dz)/(0.8*boxsigma));
                          ey=ey1 + ey2 - kbox*(1/0.8)^12-emin;
                          ey1=step(dy-0.8*boxsigma)*step(3*boxsigma-dy)*kbox*(boxsigma/dy)^12;
                          ey2=step(0.8*boxsigma-dy)*kbox*(boxsigma/(0.8*boxsigma))^12*(1 + 12*(0.8*boxsigma - dy)/(0.8*boxsigma));
                          ex=ex1 + ex2 - kbox*(1/0.8)^12-emin;
                          ex1=step(dx-0.8*boxsigma)*step(3*boxsigma-dx)*kbox*(boxsigma/dx)^12;
                          ex2=step(0.8*boxsigma-dx)*kbox*(boxsigma/(0.8*boxsigma))^12*(1 + 12*(0.8*boxsigma - dx)/(0.8*boxsigma));
                          dx=min(abs(x-0),abs(xbox-x));
                          dy=min(abs(y-0),abs(ybox-y));
                          dz=min(abs(z-0),abs(zbox-z));
                          emin=kbox*(1/3)^12;
                      """
    fix_boundary = mm.CustomExternalForce(energy_function)
    fix_boundary.addGlobalParameter('xbox',xbox*_A_to_nm)
    fix_boundary.addGlobalParameter('ybox',ybox*_A_to_nm)
    fix_boundary.addGlobalParameter('zbox',zbox*_A_to_nm)
    fix_boundary.addGlobalParameter('kbox',kbox*_kcal_to_kj)
    fix_boundary.addGlobalParameter('boxsigma',boxsigma*_A_to_nm)
    
    for ai in atoms:
        fix_boundary.addParticle(ai.index)
    fix_boundary.setForceGroup(force_group)
    
    return fix_boundary
