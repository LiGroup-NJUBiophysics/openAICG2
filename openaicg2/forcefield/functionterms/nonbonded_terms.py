import openmm as mm
from openmm import unit
import numpy as np

def oriented_dependent_Hbond_term(odHbond_info,use_pbc=False,cutoff=2.5,force_group=1):
    """a orientation dependent Lenard-Jones type potential for hydrogen bonds

    Parameters
    ----------
    odHbond_info: pd.DataFrame
        Information for all native hydrogen bonds are expressed as a N*Nine tabel.
        N row are represent the amount of Hbond. 9 columns are index of a1(atom1), 
        a2, a3, a4, a5, a6 and epsilon, sigma, cutoff(truncation distance/sigma), respectively

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default), then pbc would not apply to force.

    cutoff: float
        cutoff = truncation distance / rest length. If cutoff is None, the potential will not be trancated.

    force_group: int
        Force group.
    
    Return
    ------
    oriented_dependent_Hbond_force: Force
        OpenMM Force object
    """
    '''
        step(r-sigma)*LJ*A+(1-delta(sigma-r))*step(sigma-r)*(LJ+(1-A)*epsilon);
                                              LJ=epsilon*(5*(sigma/r)^12-6*(sigma/r)^10);
    '''
    # energy function of Hbond
    energy_func = f"""step(r-sigma)*LJ*A+(1-delta(sigma-r))*step(sigma-r)*(LJ+(1-A)*epsilon);
                                              LJ=epsilon*(5*(sigma/r)^12-6*(sigma/r)^10);
                                              r=distance(p7,p8);
                                              A=(1/(1+10*sin(theta1)^4)+1/(1+10*sin(theta2)^4))/2;  
                                              theta1=pointangle(cx1,cy1,cz1,x2,y2,z2,x5,y5,z5);
                                              theta2=pointangle(cx2,cy2,cz2,x5,y5,z5,x2,y2,z2);
                                              cx1=px1+x2;cy1=py1+y2;cz1=pz1+z2;
                                              cx2=px2+x5;cy2=py2+y5;cz2=pz1+z5;
                                              px1=(d1y*d2z-d1z*d2y);py1=(d1z*d2x-d1x*d2z);pz1=(d1x*d2y-d1y*d2x);
                                              px2=(d3y*d4z-d3z*d4y);py2=(d3z*d4x-d3x*d4z);pz2=(d3x*d4y-d3y*d4x);
                                              d1x=x1-x2;d1y=y1-y2;d1z=z1-z2;d2x=x3-x2;d2y=y3-y2;d2z=z3-z2;
                                              d3x=x4-x5;d3y=y4-y5;d3z=z4-z5;d4x=x6-x5;d4y=y6-y5;d4z=z6-z5;"""
    if cutoff != None:
        energy_func = f"""step({cutoff}*sigma - r)*""" + energy_func

    oriented_dependent_Hbond_force = mm.CustomCompoundBondForce(8,energy_func)
    oriented_dependent_Hbond_force.addPerBondParameter("epsilon")
    oriented_dependent_Hbond_force.addPerBondParameter("sigma")
    for row in odHbond_info.itertuples():
        oriented_dependent_Hbond_force.addBond([row.a1, row.a2, row.a3, row.a4, row.a5, row.a6,row.a7,row.a8], [row.epsilon, row.sigma])
    oriented_dependent_Hbond_force.setUsesPeriodicBoundaryConditions(use_pbc)
    oriented_dependent_Hbond_force.setForceGroup(force_group)
    return oriented_dependent_Hbond_force

def go_contact_term(native_contact_info, use_pbc=False, cutoff=2.5, force_group=1):
    """a typical  Lenard-Jones(12-10) potential for native contact in go model

    Parameters
    ----------
    native_contact_info: pd.DataFrame
        Information for all native contact.

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
    
    cutoff: float
        cutoff = truncated distance. If cutoff is None, the potential will not be trancated.

    force_group: int
        Force group.
    
    Return
    ------
    go_contact_force: Force
        OpenMM Force object
    """    
    if cutoff == None:
        energy_func = f'''epsilon*(5*(sigma/r)^12 - 6*(sigma/r)^10)'''
    else:
        energy_func = f'''step({cutoff}*sigma -r)*epsilon*(5*(sigma/r)^12 - 6*(sigma/r)^10)'''
    go_con_force = mm.CustomBondForce(energy_func) 
    go_con_force.addPerBondParameter('epsilon')
    go_con_force.addPerBondParameter('sigma')
    for row in native_contact_info.itertuples():
        go_con_force.addBond(row.a1, row.a2,[row.epsilon,row.sigma])
    go_con_force.setUsesPeriodicBoundaryConditions(use_pbc)
    go_con_force.setForceGroup(force_group)
    return go_con_force

def kim_hummer_term(atom_resi_type,epsilon_map, sigma_map,exclusions,extra_exclusions=None,use_pbc=False,cutoff=4,force_group=1):
    
    """a non-specific and long-range interactoin bewteen amino acid driven by their hydrophobic, aromatic, or electrostatic character, 
    which is a modified LJ type energy function and reference from "J. Mol. Biol. (2008) 375, 1416–1433". 

    Parameters
    ----------
    atom_resi_type: list, str 
        the residue name of all atoms

    epsilon_map: array
        the strength of specific residue pair

    sigma_map: array
        the rest length of specific residue pair
    
    exclusions: list
        pairs of neighbour particles whose interactions should be omitted from force and energy calculations.
    
    extra_exclusions: list
        This represents the extra particle pairs that are neglected in the calculation of the Kim-Hummer potential.

    use_pbc: bool, optional
        Whether use periodic boundary conditions. If False (default),then pbc would not apply to force.
    
    cutoff: float
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.


    force_group: int
        Force group.
    
    Return
    ------
    kim_hummer_force: Force
        OpenMM Force object
    """    
    num_type_resi = np.shape(sigma_map)[0]
    energy_func_kh = f'''step({cutoff}*sig-r)*(lj*(1-step(eps)) + step(eps)*repel);\
                          repel= f1 + f2- delta(r-sig*2^(1/6))*abs(eps);\
                          f1 = step(sig*2^(1/6)-r)*(lj+2*abs(eps));\
                          f2 = -step(r-sig*2^(1/6))*lj;\
                          lj = 4*abs(eps)*(lj12+lj6);\
                          lj12=(sig/r)^12;lj6=-(sig/r)^6;\
                          eps=epsilon(resi_type1,resi_type2);
                          sig=sigma(resi_type1,resi_type2);'''

    kim_hummer_force = mm.CustomNonbondedForce(energy_func_kh)                                  
    kim_hummer_force.addTabulatedFunction('epsilon', mm.Discrete2DFunction(num_type_resi,num_type_resi,epsilon_map.flatten()))
    kim_hummer_force.addTabulatedFunction('sigma',mm.Discrete2DFunction(num_type_resi,num_type_resi,sigma_map.flatten()))
    kim_hummer_force.addPerParticleParameter('resi_type')
    for ai in atom_resi_type:
        kim_hummer_force.addParticle([ai])
    for i_exclu in exclusions:
        kim_hummer_force.addExclusion(i_exclu[0],i_exclu[1])
    if extra_exclusions !=None:
        for i_exclu in extra_exclusions:
            if ([i_exclu[0],i_exclu[1]] == np.array(exclusions)).all(1).any() or ([i_exclu[1],i_exclu[0]] == np.array(exclusions)).all(1).any():
                continue
            kim_hummer_force.addExclusion(i_exclu[0],i_exclu[1]) 
    if use_pbc:
        kim_hummer_force.setNonbondedMethod(kim_hummer_force.CutoffPeriodic)
        kim_hummer_force.setCutoffDistance(cutoff*np.max(sigma_map))
    elif cutoff == None:
        kim_hummer_force.setNonbondedMethod(kim_hummer_force.NoCutoff)
    else:
        kim_hummer_force.setNonbondedMethod(kim_hummer_force.CutoffNonPeriodic)
        kim_hummer_force.setCutoffDistance(cutoff*np.max(sigma_map))
    kim_hummer_force.setForceGroup(force_group)
    return kim_hummer_force

def ashbaugh_hatch(atom_resi_type, epsilon, sigma_map,lambda_map, exclusions,extra_exclusions=None,use_pbc=False,cutoff=4,force_group=1):
    
    """
    ashbaugh hatch potential that is a non-specific and long-range interactoin bewteen amino acid driven by their hydrophobic, aromatic, or electrostatic character, 
    which is a modified LJ type energy function and reference from "Protein Science, 2021, 30(7): 1371-1379.". 

    Parameters
    ----------
    atom_resi_type: list, str 
        the residue type of all atoms
    
    epsilon: array
        the strength of specific residue pair

    sigma_map: array
        the rest length of specific residue pair
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    extra_exclusions: list
        This represents the extra particle pairs that are neglected in the calculation of the potential.

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
    
    cutoff: float
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.

    force_group: int
        Force group.
    
    Return
    ------
    Ash_Hat_contacts: Force
        OpenMM Force object
    """    
    num_type_resi = np.shape(sigma_map)[0]
    lj_at_cutoff = 4*epsilon*((1/cutoff)**12 - (1/cutoff)**6)
    Ash_Hat_contacts = mm.CustomNonbondedForce(f"""energy;
               energy=atom_pair_type*AH;
               AH=(f1+f2-offset)*step(4*sigma_ah-r);
               offset=lambda_ah*{lj_at_cutoff};
               f1=(lj+(1-lambda_ah)*{epsilon})*step(2^(1/6)*sigma_ah-r);
               f2=lambda_ah*lj*step(r-2^(1/6)*sigma_ah);
               lj=4*{epsilon}*(lj12 + lj6);
               lj12=(sigma_ah/r)^12;
               lj6=-(sigma_ah/r)^6;
               sigma_ah=sigma_ah_map(resi_type1, resi_type2);
               lambda_ah=lambda_ah_map(resi_type1, resi_type2);
               atom_pair_type=atom_type_map(atom_type1,atom_type2);
               """)
    discrete_2d_sigma_ah_map = mm.Discrete2DFunction(num_type_resi, num_type_resi, sigma_map.ravel().tolist())
    discrete_2d_lambda_ah_map = mm.Discrete2DFunction(num_type_resi, num_type_resi, lambda_map.ravel().tolist())
    Ash_Hat_contacts.addTabulatedFunction('sigma_ah_map', discrete_2d_sigma_ah_map)
    Ash_Hat_contacts.addTabulatedFunction('lambda_ah_map', discrete_2d_lambda_ah_map)
    Ash_Hat_contacts.addPerParticleParameter('resi_type')
    for ai in atom_resi_type:
        Ash_Hat_contacts.addParticle([ai])
    for i_exclu in exclusions:
        Ash_Hat_contacts.addExclusion(i_exclu[0],i_exclu[1])
    if extra_exclusions != None:
        for i_exclu in extra_exclusions:
            if ([i_exclu[0],i_exclu[1]] == np.array(exclusions)).all(1).any() or ([i_exclu[1],i_exclu[0]] == np.array(exclusions)).all(1).any():
                continue
            Ash_Hat_contacts.addExclusion(i_exclu[0],i_exclu[1]) 
    if use_pbc:
        Ash_Hat_contacts.setNonbondedMethod(Ash_Hat_contacts.CutoffPeriodic)
        Ash_Hat_contacts.setCutoffDistance(cutoff*np.max(sigma_map))
    elif cutoff == None:
        Ash_Hat_contacts.setNonbondedMethod(Ash_Hat_contacts.NoCutoff)
    else:
        Ash_Hat_contacts.setNonbondedMethod(Ash_Hat_contacts.CutoffNonPeriodic)
        Ash_Hat_contacts.setCutoffDistance(cutoff*np.max(sigma_map))
    Ash_Hat_contacts.setForceGroup(force_group)
    return Ash_Hat_contacts

def excluded_term(atom_types,epsilon_map,sigma_map,exclusions,extra_exclusions=None,use_pbc=False,cutoff=2.5,force_group=1):
    """
    An excluded interaction between atoms is implemented to prevent overlap.

    Parameters
    ----------
    atom_type: list, int 
        Atom type.
    
    epsilon_map: array
        the strength of specific residue pair

    sigma_map: array
        the rest length of specific residue pair
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    extra_exclusions: list
        This represents the extra particle pairs that are neglected in the calculation of the potential.

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.

    cutoff: float
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.
    
    force_group: int
        Force group.
    
    Return
    ------
    excluded_force: Force
        OpenMM Force object    
    """
    num_type_atom = np.shape(epsilon_map)[0]
    excluded_force = mm.CustomNonbondedForce(f'''step({cutoff}*sig-r)*eps*(sig/r)^12;
                                              eps=epsilon(atom_type1,atom_type2);
                                              sig=sigma(atom_type1,atom_type2);''')#;ex_resi=(step(abs(ex_resi1-ex_resi2)-3))') 
    excluded_force.addTabulatedFunction('epsilon', mm.Discrete2DFunction(num_type_atom,num_type_atom,epsilon_map.flatten()))
    excluded_force.addTabulatedFunction('sigma', mm.Discrete2DFunction(num_type_atom,num_type_atom,sigma_map.flatten()))
    excluded_force.addPerParticleParameter('atom_type')
    for ai in atom_types:
        excluded_force.addParticle([ai])
    for i_exclu in exclusions:
        excluded_force.addExclusion(i_exclu[0],i_exclu[1])
    if extra_exclusions != None:
        for i_exclu in extra_exclusions:
            if ([i_exclu[0],i_exclu[1]] == np.array(exclusions)).all(1).any() or ([i_exclu[1],i_exclu[0]] == np.array(exclusions)).all(1).any():
                continue
            excluded_force.addExclusion(i_exclu[0],i_exclu[1]) 
    if use_pbc:
        excluded_force.setNonbondedMethod(excluded_force.CutoffPeriodic)
        excluded_force.setCutoffDistance(cutoff*np.max(sigma_map))
    elif cutoff == None:
        excluded_force.setNonbondedMethod(excluded_force.NoCutoff)
    else:
        excluded_force.setNonbondedMethod(excluded_force.CutoffNonPeriodic)
        excluded_force.setCutoffDistance(cutoff*np.max(sigma_map))
    excluded_force.setForceGroup(force_group)
    return excluded_force

def kh_and_ex_term(resi_type_list,atom_type_list,epsilon_map,sigma_map,atom_type_map,exclusions,extra_exclusions=None,use_pbc=False,cutoff_kh=2.5,cutoff_ex=2.0,force_group=1):
    """
    This function implements kim-hummer and excluded potential.

    Parameters
    ----------
    resi_type_list: int
        residual type.

    atom_type_list: int 
        Atom type.
    
    epsilon_map: array
        the strength of specific atom pair

    sigma_map: array
        the rest length of specific atom pair
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    extra_exclusions: list
        This represents the extra particle pairs that are neglected in the calculation of the potential.

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.

    cutoff_kh: float
        cutoff for kim-hummer potential.
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.
    
    cutoff_ex: float
        cutoff for excluded potential.
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.
 
    force_group: int
        Force group.
    
    Return
    ------
    KH_ex_force: Force
        OpenMM Force object   
    """
    num_type_resi = np.shape(epsilon_map)[0] - 1
    num_atom_type = np.shape(atom_type_map)[0]
    energy_func_part1 = """  rep= f1 + f2-delta(r-sig*2^(1/6))*abs(eps);
                                   f1 = step(sig*2^(1/6)-r)*(lj+2*abs(eps));
                                   f2 = -step(r-sig*2^(1/6))*lj;
                                   lj=4*abs(eps)*(lj12+lj6);
                                   ex=abs(eps)*lj12; 
                                   lj12=(sig/r)^12;lj6=-(sig/r)^6;
                                   eps=epsilon(resi_type1,resi_type2);
                                   sig=sigma(resi_type1,resi_type2);
                                   atom_type=atom_type_map(atom_type1,atom_type2);"""


    if cutoff_kh != None and cutoff_ex!=None:
        KH_and_ex_ener_func =  f"""(1-atom_type)*step({cutoff_kh}*sig-r)*(lj*(1-step(eps)) + step(eps)*rep) + atom_type*step({cutoff_ex}*sig-r)*ex;""" + energy_func_part1
    elif cutoff_kh == None or cutoff_ex == None:
        KH_and_ex_ener_func = """(1-atom_type)*(lj*(1-step(eps)) + step(eps)*rep) + atom_type*ex;""" + energy_func_part1


    KH_ex_force = mm.CustomNonbondedForce(KH_and_ex_ener_func)
    KH_ex_force.addTabulatedFunction('epsilon', mm.Discrete2DFunction(int(num_type_resi+1),int(num_type_resi+1),epsilon_map.flatten()))
    KH_ex_force.addTabulatedFunction('sigma',mm.Discrete2DFunction(int(num_type_resi+1),int(num_type_resi+1),sigma_map.flatten()))
    KH_ex_force.addTabulatedFunction('atom_type_map',mm.Discrete2DFunction(num_atom_type,num_atom_type,atom_type_map.flatten()))
    KH_ex_force.addPerParticleParameter('resi_type')
    KH_ex_force.addPerParticleParameter('atom_type')
    for i,resi_type in enumerate(resi_type_list):
        KH_ex_force.addParticle([resi_type,atom_type_list[i]])
    for i_exclu in exclusions:
        KH_ex_force.addExclusion(i_exclu[0],i_exclu[1])
    if extra_exclusions !=None:
        for i_exclu in extra_exclusions:
            if ([i_exclu[0],i_exclu[1]] == np.array(exclusions)).all(1).any() or ([i_exclu[1],i_exclu[0]] == np.array(exclusions)).all(1).any():
                continue
            KH_ex_force.addExclusion(i_exclu[0],i_exclu[1]) 
    if use_pbc:
        KH_ex_force.setNonbondedMethod(KH_ex_force.CutoffPeriodic)
        KH_ex_force.setCutoffDistance(cutoff_kh*np.max(sigma_map))
    elif cutoff_kh ==None or cutoff_ex == None:
        KH_ex_force.setNonbondedMethod(KH_ex_force.NoCutoff)
    else:
        KH_ex_force.setNonbondedMethod(KH_ex_force.CutoffNonPeriodic)
        KH_ex_force.setCutoffDistance(cutoff_kh*np.max(sigma_map))
    KH_ex_force.setForceGroup(force_group)
    return KH_ex_force

def debye_Huckel_bond_form(charge,exclusions,dieletric_constant=80,ion_strength=0.02,T=300,use_pbc=False,cutoff=20,force_group=1):
    '''a Debye-Huckel potential is acheived by CustomBondForce
    
    Parameters
    ----------
    charge: pd.DataFrame
        include the index and charge quatity of pair atom 
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    dieletric_constant: float
        dieletric constant

    ion_strength: float, unit is mol/L
        salt-conentratoin,the default value 0.02 is correspond 20mM/L NaCl and ion_strength I=1 for 1M/L NaCl
    
    T: float, unit.kelvin
        temperatur
    
    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.

    cutoff: float
        ele_cutoff = (truncated distane/lamdaD). lamdaD is so-called Debye length

    force_group: int
        Force group.
        
    '''

    epsilon_0 = 8.8542e-12 #unit.farad/unit.meter c/(m*V); electric contant
    lamdaD = (epsilon_0*1.380649e-23*1e17*dieletric_constant*T/(2*6.02e23 *ion_strength* (1.602e-19)**2))**0.5/10 #*unit.nanometer # Debye length
    coeff_eps = (1.602e-19)**2*1e9/(4*np.pi*epsilon_0*dieletric_constant*6.9478e-21)
    energy_func = '''coeff_1*charge1*charge2*exp(-r/lamdaD)/r'''
    if cutoff != None:
        energy_func = f'''step({lamdaD*cutoff}-r)*''' + energy_func 
    dh_force = mm.CustomBondForce(energy_func)
    dh_force.addPerBondParameter('charge1')
    dh_force.addPerBondParameter('charge2')
    dh_force.addGlobalParameter('coeff_1',coeff_eps*unit.kilocalorie_per_mole*unit.nanometer)
    dh_force.addGlobalParameter('lamdaD',lamdaD*unit.nanometer)
    for i_ele in charge:
        # set exclusion for exluding neighbor atom
        if ([i_ele[0],i_ele[1]] == np.array(exclusions)).all(1).any() or ([i_ele[1],i_ele[0]] == np.array(exclusions)).all(1).any():
           continue
        dh_force.addBond(i_ele[0],i_ele[1],i_ele[2:])
    dh_force.setUsesPeriodicBoundaryConditions(use_pbc)
    dh_force.setForceGroup(force_group)
    return dh_force

def debye_Huckel_nonbonded_form(charge,exclusions,dieletric_constant=80,ion_strength=0.02,T=300,use_pbc=False,cutoff=20,force_group=1):
    '''a Debye-Huckel potential is acheived by CustomNonbondedForce
    
    Parameters
    ----------
    charge: list
        charge of each atoms
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    dieletric_constant: float
        dieletric constant

    ion_strength: float, unit is mol/L
        salt-conentratoin,the default value 0.02 is correspond 20mM/L NaCl and ion_strength I=1 for 1M/L NaCl
    
    T: float, unit.kelvin
        temperature

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
        
    cutoff: float
        ele_cutoff = (truncated distane/lamdaD). lamdaD is so-called Debye length
    
    force_group: int
        Force group.
    
    Return
    ------
    deb_huc_force: Force
        OpenMM Force object  
    '''

    epsilon_0 = 8.8542e-12 #unit.farad/unit.meter c/(m*V); electric contant
    lamdaD = (epsilon_0*1.380649e-23*1e17*dieletric_constant*T/(2*6.02e23 *ion_strength* (1.602e-19)**2))**0.5/10 #*unit.nanometer # Debye length
    coeff_eps = (1.602e-19)**2*1e9/(4*np.pi*epsilon_0*dieletric_constant*6.9478e-21) 
    dh_force_nonbonded = mm.CustomNonbondedForce('step(ele_cutoff-r)*coeff_1*charge1*charge2*exp(-r/lamdaD)/r')
    dh_force_nonbonded.addPerParticleParameter('charge1')
    dh_force_nonbonded.addPerParticleParameter('charge2')
    dh_force_nonbonded.addGlobalParameter('coeff_1',coeff_eps*unit.kilocalorie_per_mole*unit.nanometer)
    dh_force_nonbonded.addGlobalParameter('lamdaD',lamdaD*unit.nanometer)
    for iq in charge:
        dh_force_nonbonded.addParticle(iq)
    for i_exclu in exclusions:
        dh_force_nonbonded.addExclusion(i_exclu[0],i_exclu[1])
    if use_pbc:
        dh_force_nonbonded.setNonbondedMethod(dh_force_nonbonded.CutoffPeriodic)
        dh_force_nonbonded.setCutoffDistance(lamdaD*ele_cutoff*unit.nanometer)
    elif cutoff == None:
        dh_force_nonbonded.setNonbondedMethod(dh_force_nonbonded.NonCutoff)
    else:
        dh_force_nonbonded.setNonbondedMethod(dh_force_nonbonded.CutoffNonPeriodic)
        dh_force_nonbonded.setCutoffDistance(lamdaD*cutoff*unit.nanometer)
    dh_force_nonbonded.setForceGroup(force_group)
    return dh_force_nonbonded
