import numpy as np

_kcal_to_kj = 4.1840

class parser_flp_para(object):

    def __init__(self):
        """
        """
        self.FBA_MIN_ANG = 1.31
        self.FBA_MAX_ANG = 2.87
        self.FBA_MIN_ANG_FORCE = -30.0# * _kcal_to_kj
        self.FBA_MAX_ANG_FORCE = 30.0# * _kcal_to_kj
        self.flp_bond_dihd_params = {}
        self.bond_ang_x = []
        self.bond_ang_y = {}
        self.bond_ang_y2 = {}
        self.flp_bond_ang_params = {}

    def read_flexible_local_para(self,path):
        '''
        read the parameter file of flexible local potential

        Parameters
        ----------
        path : str, 
            The path of flexible local potential parameter file.ls
        
        Return
        ------
        flp_dihd: dictionary, float
            This parameters are used to describe the virtural dihedral angle in the disorder region.
            Each residue pair has a specific set of parameters. From the first column onwards, These 
            parameters correspond to the constant term C and coefficients c1, s1, c2, s2, c3, and s3 in 
            the Fourier fitting formula. Fourier fitting formula: C+c1*cos(theta)+s1*sin(theta)+ c2*cos(2*theta)+s2*sin(2*theta)+
            c3*cos(3*theta)+s3*sin(3*theta). 
        
        flp_ang_x: list, float
            the angle interval of cubic spline

        flp_ang_y: dictionary, float, kj/(mol*radin)
            The y-value at the endpoint of the interval depends on the residue.

        flp_ang_y2: dictionary, float, kj/(mol*radin^2)
            The second derivative of spline interpolation corresponds to the angle at the interval endpoints.
        
        '''
        with open(path,'r') as read_flp:
            for line in read_flp:
                line = line.strip('\n')
                if line == '>>>>':
                    to_read = 0
                    continue
                if len(line) == 0:
                    continue
                if line == '<<<< dihedral_angle':
                    to_read = 'dihd'
                elif to_read == 'dihd':
                    line = line.split()
                    self.flp_bond_dihd_params[line[0] + line[1]] = np.array(line[2:]).astype(np.float64) #* _kcal_to_kj
                elif line == '<<<< bond_angle_x':
                    to_read = 'angle x'
                elif to_read == 'angle x':
                    line = line.split()
                    self.bond_ang_x.append(float(line[1]))
                elif line == '<<<< bond_angle_y':
                    to_read = 'angle y'
                elif to_read == 'angle y':
                    line = line.split()
                    self.bond_ang_y[line[0]] = np.array(line[1:]).astype(np.float64) #* _kcal_to_kj
                elif line == '<<<< bond_angle_y2':
                    to_read = 'angle y2'
                elif to_read == 'angle y2':
                    line = line.split()
                    self.bond_ang_y2[line[0]] = np.array(line[1:]).astype(np.float64)# * _kcal_to_kj
    
    # Get parameter for correcting the energy of flexible local angle 
    def cubic_spline(self,theta,theta_lo,theta_hi,y_lo,y_hi,y2_lo,y2_hi):
        """
        The Cubic spline function.

        Parameters
        ----------
        theta: float
            angle

        theta_lo: float
            the low limit of theta

        theta_hi: float
            the high limit of theta
        
        y_lo: float
            the energy value correspond to theta_lo

        y_hi: float
            the energy value correspond to theta_hi
        
        y2_lo: float
            the second-order derivative value of energy correspond to theta_lo
        
        y2_hi: float
            the second-order derivative value of energy correspond to theta_hi
        
        return
        ------
        energy: float
            the energy of flexible local angle
        """
        lk = theta_hi - theta_lo
        a = (theta_hi - theta) / lk
        b = (theta - theta_lo) / lk    
        return (a**3 - a)*lk**2*y2_lo/6 + (b**3-b)*lk**2*y2_hi/6 + b * y_hi + a * y_lo

    def diff_cubic_spline(self,theta,theta_lo,theta_hi,y_lo,y_hi,y2_lo,y2_hi):
        lk = theta_hi - theta_lo
        a = (theta_hi - theta) / lk
        b = (theta - theta_lo) / lk    
        return (1-3*a**2)*lk*y2_lo/6 + (3*b**2-1)*lk*y2_hi/6 + y_hi/lk - y_lo/lk

    def correct_flex_ang_force_para(self,num_bin: int,para,stepsize=1e-4):
        bond_angle_x = para[:num_bin+1]
        bond_angle_y = para[num_bin+1:2*(num_bin+1)]
        bond_angle_y2 = para[2*(num_bin+1):]
        # initial the boudary parameter
        min_theta = self.FBA_MIN_ANG
        max_theta = self.FBA_MIN_ANG
        centre_theta = (self.FBA_MAX_ANG - self.FBA_MIN_ANG)/2
        min_theta_energy = self.cubic_spline(min_theta,bond_angle_x[0],bond_angle_x[1],bond_angle_y[0],bond_angle_y[1],bond_angle_y2[0],bond_angle_y2[1])
        for i in range(num_bin):
            theta = np.arange(bond_angle_x[i],bond_angle_x[i+1],stepsize)
            energy = self.cubic_spline(theta,bond_angle_x[i],bond_angle_x[i+1],bond_angle_y[i],bond_angle_y[i+1],bond_angle_y2[i],bond_angle_y2[i+1])
            force = self.diff_cubic_spline(theta,bond_angle_x[i],bond_angle_x[i+1],bond_angle_y[i],bond_angle_y[i+1],bond_angle_y2[i],bond_angle_y2[i+1])
            if i == 0:
                minimum_energy = np.min(energy)
            else:
                mini_ener_i = np.min(energy)
                minimum_energy = np.min([minimum_energy,mini_ener_i])
            index_mini  = np.argwhere(force < self.FBA_MIN_ANG_FORCE)
            if len(index_mini) != 0:
               min_theta = theta[index_mini[-1]]
               min_theta_energy = energy[index_mini[-1]]
            index_max = np.argwhere(force > self.FBA_MAX_ANG_FORCE)
            if len(index_max)  != 0 and max_theta == self.FBA_MIN_ANG:
                max_theta = theta[index_max[0]]
                if max_theta > centre_theta:
                   max_theta_energy = energy[index_max[0]]
        bond_angle_x = np.array(bond_angle_x)
        bond_angle_y = np.array(bond_angle_y) * _kcal_to_kj
        bond_angle_y2 = np.array(bond_angle_y2) * _kcal_to_kj
        boundary_para = np.array([min_theta,min_theta_energy*_kcal_to_kj,max_theta,max_theta_energy*_kcal_to_kj,minimum_energy*_kcal_to_kj],dtype=object)
        corr_para = np.concatenate((bond_angle_x,bond_angle_y,bond_angle_y2,boundary_para),axis=0)
        return corr_para
    # set the flexible local energy corr? for dihedral interaction of backbond
    def  flexi_dihd_energy(self,theta,para):
        return  para[0] + para[1]*np.cos(theta) + para[2]*np.sin(theta) + para[3]*np.cos(2*theta) + para[4]*np.sin(2*theta) \
                + para[5]*np.cos(3*theta) + para[6]*np.sin(3*theta)
    
    def set_flex_dihd_corr(self,para,stepsize=1e-4):
        theta_lo = -np.pi
        theta_hi = np.pi
        theta = np.arange(theta_lo,theta_hi,stepsize)
        ener = self.flexi_dihd_energy(theta,para)
        para_all = np.hstack((para,np.min(ener)))
        return para_all*_kcal_to_kj
        
    def get_corr_flex_ang_para(self,path):
        """
        Obtain the corrected flexible local potential parameters for the bond angle and dihedral angle of any residue
        """
        self.read_flexible_local_para(path)
        residue_list = self.bond_ang_y.keys()
        for i_resi in residue_list:
            num_interval = len(self.bond_ang_x) - 1
            flp_para = np.concatenate((self.bond_ang_x,self.bond_ang_y[i_resi],self.bond_ang_y2[i_resi]),axis=0)
            self.flp_bond_ang_params[i_resi] = self.correct_flex_ang_force_para(num_interval,flp_para)
        residue_pair = self.flp_bond_dihd_params.keys()
        for i_resi_pair in residue_pair:
            resi_pair_para = self.set_flex_dihd_corr(self.flp_bond_dihd_params[i_resi_pair])
            self.flp_bond_dihd_params[i_resi_pair] = resi_pair_para


            


