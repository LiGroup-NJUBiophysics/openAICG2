import numpy as np 
from  openmm import unit

_kcal_to_kj = 4.1840
_A_to_nm = 0.1
def parser_exv_params(path):
    """
    Reading excluded potential parameters.

    Parameters
    ---------
    path: str
       The name of path
    
    return
    ------
    exv_params: dictionary
        exv_params include the following keywords: "SIGMA", "xv_cutoff" and "exv_coef".
    """

    exv_params = {}
    with open(path,'r') as read_exv:
        for line in read_exv:
            line.strip()
            if 'SIGMA' in line:
                line = line.split()
                exv_params[line[1]] = float(line[2]) * _A_to_nm
            elif 'exv_cutoff' in line:
                line = line.split()
                exv_params[line[0]] = float(line[2]) * _A_to_nm
            elif 'exv_coef' in line:
                line = line.split()
                exv_params[line[0]] = float(line[2]) * _kcal_to_kj
    return exv_params
                
                
