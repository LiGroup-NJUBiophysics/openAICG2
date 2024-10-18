import pandas as pd
import openmm as mm
from openmm import unit

CUTOFF_UNDER_EXP = -50
def harmonic_angle_term(angle_info,use_pbc=False,force_group=1):
    """
    Harmonic angle force term

    Parameters
    ----------
    angle_info : pd.DataFrame
        Information about all Hamonic angles

    use_pbc : bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to Harmonic angle force. 
    
    force_group : int
        Force group
    
    Return
    ------
    angle_force : Force
        OpenMM harmonic angle force object 
    """
    angle_force = mm.HarmonicAngleForce()
    for row in angle_info.itertuples():
        angle_force.addAngle(row.a1, row.a2, row.a3, row.natang,row.k)
    angle_force.setUsesPeriodicBoundaryConditions(use_pbc)
    angle_force.setForceGroup(force_group)
    return angle_force

def aicg13_angle_term(aicg13_ang_info,use_pbc=False,force_group=1):
    """
    aicg13 interaction is structure-based local contact potential to describe specific local interactions of the given protein structure. 
    
    Parameters
    ----------
    aicg13_info : pb.DataFrame
        Information for all aicg13 angle interaction.
    
    use_pbc : bool, optional 
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to aicg13 interaction.
    
    force_group : int
        Force group
    
    Return
    ------
    aicg13_ang_force : Force
        Openmm force object
    """
    #CUTOFF_UNDER_EXP = -50
    #*step(gex-{CUTOFF_UNDER_EXP})
    aicg13_ang_force = mm.CustomCompoundBondForce(3,f"-epsilon13 * exp(gex);gex=-(r13-r13_0)^2/(2*width13^2);r13=distance(p3,p1)")
    aicg13_ang_force.addPerBondParameter("epsilon13")
    aicg13_ang_force.addPerBondParameter("r13_0")
    aicg13_ang_force.addPerBondParameter("width13")
    for row in aicg13_ang_info.itertuples():
        aicg13_ang_force.addBond([row.a1, row.a2, row.a3], [row.epsilon, row.r0, row.width])
    aicg13_ang_force.setUsesPeriodicBoundaryConditions(use_pbc)
    aicg13_ang_force.setForceGroup(force_group)
    return aicg13_ang_force 

def flex_angle_term(flex_ang_info,use_pbc=False,force_group=1):
    """
     The flexible local potential was constructed by analyzing loop structures in protein structure database to enchance angle flexibility.

    Parameters
    ----------
    flex_ang_info : pb.DataFrame
        Information of flexible local potential for the virtual bond angles.
    
    use_pbc : bool, optional 
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.
    
    force_group : int
        Force group.
    
    Return
    ------
    flex_ang_force : Force
        Openmm force object.
    """
    
    FBA_MIN_ANG_FORCE = -30.0
    FBA_MAX_ANG_FORCE = 30.0
    num_interval = 9#len(flex_angle_x) - 1
    theta_samp_x_name = ['theta_%d'%i for i in range(num_interval+1)]
    theta_y2_name = ['y2_%d'%i for i in range(num_interval+1)]
    theta_y_name = ['y_%d'%i for i in range(num_interval+1)]
    #    a = (theta_hi - theta) / lk
    #    b = (theta - theta_lo) / lk
    #    return (a**3 - a)*lk**2*y2_lo/6 + (b**3-b)*lk**2*y2_hi/6 + b * y_hi + a * y_lo
    for ind in range(num_interval):
        equa_fisrt_part = "(((%s - theta)^3/(%s-%s) - (%s - theta)*(%s-%s))*%s/6 + " \
                          %(theta_samp_x_name[ind+1],theta_samp_x_name[ind+1],
                            theta_samp_x_name[ind],theta_samp_x_name[ind+1],
                            theta_samp_x_name[ind+1],theta_samp_x_name[ind],theta_y2_name[ind])
        equa_sec_part = "((theta - %s)^3/(%s-%s) - (theta- %s)*(%s-%s))*%s/6 + "\
                        %(theta_samp_x_name[ind],theta_samp_x_name[ind+1],
                          theta_samp_x_name[ind],theta_samp_x_name[ind],
                          theta_samp_x_name[ind+1],theta_samp_x_name[ind],theta_y2_name[ind+1])
        equa_thr_part = "(theta - %s)*%s/(%s-%s) + (%s - theta) * %s/(%s-%s))"\
                        %(theta_samp_x_name[ind],theta_y_name[ind+1],
                          theta_samp_x_name[ind+1],theta_samp_x_name[ind],
                          theta_samp_x_name[ind+1],theta_y_name[ind],
                          theta_samp_x_name[ind+1],theta_samp_x_name[ind])
        equa_interv_part = "*step(theta-%s)*step(%s-theta)"%(theta_samp_x_name[ind],theta_samp_x_name[ind+1])
        if ind == 0:
           ener_expression = equa_fisrt_part + equa_sec_part + equa_thr_part + equa_interv_part
        else:
           rid_boundary_val = "-delta(theta-%s)*%s"%(theta_samp_x_name[ind],theta_y_name[ind])
           interval_expre_i = '+' + equa_fisrt_part + equa_sec_part + equa_thr_part +  equa_interv_part + rid_boundary_val
           ener_expression += interval_expre_i   

    FBA_MIN_ANG_FORCE_UNIT = FBA_MIN_ANG_FORCE * 4.1840#* unit.kilocalorie_per_mole/unit.radian
    FBA_MAX_ANG_FORCE_UNIT = FBA_MAX_ANG_FORCE * 4.1840#* unit.kilocalorie_per_mole/unit.radian
    ener_expression = 'step(theta - fba_min_th)*step(fba_max_th-theta)*('+ener_expression + ')' +\
                      f'+step(fba_min_th-theta) * ({FBA_MIN_ANG_FORCE_UNIT}*theta + fba_min_th_ener-{FBA_MIN_ANG_FORCE_UNIT}* fba_min_th)' +\
                      f'+ step(theta-fba_max_th) * ({FBA_MAX_ANG_FORCE_UNIT}*theta + fba_max_th_ener - {FBA_MAX_ANG_FORCE_UNIT} * fba_max_th)' + \
                       '- delta(theta-fba_min_th) * fba_min_th_ener - delta(fba_max_th - theta) * fba_max_th_ener - corr_tot_fba_ener'
    # create a Custom angle force for flexible angle potential   
    bound_ener_corr_para_name = ['fba_min_th','fba_min_th_ener','fba_max_th','fba_max_th_ener','corr_tot_fba_ener']
    flex_angle_force = mm.CustomAngleForce(ener_expression)
    for i in range(len(theta_samp_x_name)):
        flex_angle_force.addPerAngleParameter(theta_samp_x_name[i])
    for i in range(len(theta_y_name)):
        flex_angle_force.addPerAngleParameter(theta_y_name[i])
    for i in range(len(theta_y2_name)):
        flex_angle_force.addPerAngleParameter(theta_y2_name[i])
    for i in bound_ener_corr_para_name:
        flex_angle_force.addPerAngleParameter(i)
    for _, row in flex_ang_info.iterrows():
        a1 = int(row['a1'])
        a2 = int(row['a2'])
        a3 = int(row['a3'])
        para = row[3:]
        flex_angle_force.addAngle(a1, a2, a3,para)
    flex_angle_force.setUsesPeriodicBoundaryConditions(use_pbc)
    flex_angle_force.setForceGroup(force_group)
    return flex_angle_force
           
