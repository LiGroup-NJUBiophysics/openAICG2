import numpy as np
from openmm import unit


def parser_kh_params(path,kh_model_symbol,kh_epsilon_scale,T):
    """
    Reading kim-hummer potentail parameter file.

    Parameters
    ----------
    path: str
        Path name of parameter file.

    kh_model_symbol: str
        The Kim-Hummer potential has six model parameters, which you can choose from A, B, C, D, E, and F.

    kh_epsilon_scale: float
        The parameter scales the strength of lambda to generate epsilon value.

    T: float
        Current simulation temperature.

    Return
    ------
    resi_type: string list
        20 residue type
    
    epsilon_KH: numpy.array 
        The parameter array consists of a 21x21 matrix for epsilon values. The 20x20 section of the array 
        represents the parameter matrix for the side chain bead (CB) residue types, while the 21st row and 
        21st column represent the parameters for the central atom (CA).
    
    sigma_KH: numpy.array
        The parameter array consists of a 21x21 matrix for sigma values. The 20x20 section of the array 
        represents the parameter matrix for the side chain bead (CB) residue types, while the 21st row and 
        21st column represent the parameters for the central atom (CA).

    """
    _kcal_to_kj = 4.184
    _A_to_nm = 0.1
    kBT = unit.BOLTZMANN_CONSTANT_kB * T *unit.kelvin/(6.9478e-21*unit.joule) # unit: kcal/mol
    # load and arrange parameter array for Kim-Hummer potential
    kh_para_all = {}
    with open(path,'r') as read_kh:
        for line in read_kh:
            line = line.strip()
            if line == "<<<< kh_para" or line == "<<<< lamda_e0" or line == "<<<< sigma_kh":
                line = line.split()
                kh_para_all[line[1]] = []
                para_key = line[1]
            elif line == '>>>>':
                continue
            else:
                line = line.split()
                kh_para_all[para_key].append(line)
    eij_kh = {}
    for kh_i in kh_para_all['kh_para']:
        eij_kh[kh_i[0]+kh_i[1]] = float(kh_i[2])
    kh_model_para = np.array(kh_para_all['lamda_e0']).astype(np.float64)
    model = dict()
    model_symbol = ['A','B','C','D','E','F']
    for im,mparams in enumerate(kh_model_para):
        model[model_symbol[im]] = mparams[1:]
    lambda_e0 = model[kh_model_symbol]
    # construct sigma and epsilon TabulatedFunction array for exclude force and kim hummer force
    name_resi = [sig[0] for sig in kh_para_all['sigma_kh']]#list(sigma_KH.keys())
    num_type_resi = len(name_resi)
    sigma_value = [float(sig[1]) for sig in kh_para_all['sigma_kh']]
    epsilon_KH = np.zeros((num_type_resi+1,num_type_resi+1),dtype=float)
    sigma_KH = np.zeros((num_type_resi+1,num_type_resi+1),dtype=float)
    resi_type = {}
    for i in range(num_type_resi):
        resi_type[name_resi[i]] = i
        for j in range(num_type_resi):
            keys_res = name_resi[i] + name_resi[j]
            if keys_res in eij_kh.keys():
                pass
            else:
                keys_res = name_resi[j] + name_resi[i]
            epsilon_ij = lambda_e0[0] * (eij_kh[keys_res]-lambda_e0[1])*kBT
            if abs(epsilon_ij) < 0.001:
                if epsilon_ij < 0:
                    epsilon_ij = -0.001
                elif epsilon_ij >= 0:
                    epsilon_ij = 0.001
            epsilon_KH[i,j] = kh_epsilon_scale * epsilon_ij * _kcal_to_kj
            sigma_KH[i,j] = _A_to_nm*(sigma_value[i] + sigma_value[j])/2
    return resi_type, epsilon_KH, sigma_KH
