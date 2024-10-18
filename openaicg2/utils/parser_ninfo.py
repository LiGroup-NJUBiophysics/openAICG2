import numpy as np
import pandas as pd

_A_to_nm = 0.1
_kcal_to_kj = 4.1840
class ParserNinfo():
     def __init__(self):
        """
        Initialize
        """
        self.forcefield_term = ['protein_bonds', 'protein_harmonic_angles','protein_aicg13_angles',
                                'protein_native_dihd', 'protein_aicg_dihd', 'protein_native_pair']

     def convert_str_to_number_array(self,para_array,int_ini_idx,int_end_idx,float_ini_idx,float_end_idx):
          para_array = np.array(para_array)
          idx_array = para_array[:,int_ini_idx:int_end_idx].astype(np.int64) - int(1) 
          kpara_array = para_array[:,float_ini_idx:float_end_idx].astype(np.float64)
          para_array = [idx_array,kpara_array]
          return para_array

     def parser_ninfo(self,ninfo_file_path):
          """
          The method parser native information file.

          Parameter: str
             the path of native information file.
          """

          force_field_para = {}
          with open(ninfo_file_path,'r') as read_f:
               info = []
               for line in read_f:
                    line = line.strip('\n')
                    if line == '>>>>':
                         force_field_para[keys] = info
                         info = []
                    line = line.split()
                    if len(line) == 0:
                         continue
                    elif line[0] == 'bond' and int(line[2]) == 1:
                         keys = line[0]
                         info.append(line)
                    elif line[0] == 'angl' and int(line[2]) == 1:
                         keys = line[0]
                         info.append(line)
                    elif line[0] == 'aicg13' and int(line[2]) == 1:
                         keys = line[0]
                         info.append(line)
                    elif line[0] == 'dihd' and int(line[2]) == 1:
                         keys = line[0]
                         info.append(line)
                    elif line[0] == 'aicgdih' and int(line[2]) == 1:
                         keys = line[0]
                         info.append(line)
                    elif line[0] == 'contact' and int(line[2]) == 1 and int(line[3]) == 1:
                         keys = line[0]
                         info.append(line)
                    elif line[0] == 'contact' and int(line[2]) == 1 and int(line[3]) == 2:
                         keys = line[0]
                         info.append(line)
          if 'bond' in force_field_para:
             self.protein_bonds = self.convert_str_to_number_array(force_field_para['bond'],1,8,8,12)
          if 'angl' in force_field_para:
             self.protein_harmonic_angles = self.convert_str_to_number_array(force_field_para['angl'],1,10,10,14)
          if 'aicg13' in force_field_para:
             self.protein_aicg13_angles = self.convert_str_to_number_array(force_field_para['aicg13'],1,10,10,15)
          if 'dihd' in force_field_para:
             self.protein_native_dihd = self.convert_str_to_number_array(force_field_para['dihd'],1,12,12,17)
          if 'aicgdih' in force_field_para:
             self.protein_aicg_dihd = self.convert_str_to_number_array(force_field_para['aicgdih'],1,12,12,17)
          if 'contact' in force_field_para:
             self.protein_native_pair = self.convert_str_to_number_array(force_field_para['contact'],1,8,8,12)
     
     def bonds_array_to_pd(self):
          """
          To make protein bonds array to tabular format

          """
          if hasattr(self,'protein_bonds'):
               idx_bonds = self.protein_bonds[0][:,5:7]
               bd_nat = self.protein_bonds[1][:,0] * _A_to_nm
               coef_bd = self.protein_bonds[1][:,3] * _kcal_to_kj * 100 * 2
               pd_idx_bonds = pd.DataFrame(idx_bonds, columns=['a1','a2'])
               pd_bd_nat = pd.DataFrame(bd_nat, columns=['r0'])
               pd_coef_bd = pd.DataFrame(coef_bd, columns=['k'])
               self.protein_bonds = pd.concat([pd_idx_bonds,pd_bd_nat,pd_coef_bd],axis=1)
     
     def harm_ang_array_to_pd(self):
          """
          To make harmonic angles array to tabular format

          """
          if hasattr(self,'protein_harmonic_angles'):
               idx_harm_ang = self.protein_harmonic_angles[0][:,6:9]
               nat_ang = self.protein_harmonic_angles[1][:,0] * np.pi / 180
               coef_ang = self.protein_harmonic_angles[1][:,3] * _kcal_to_kj * 2
               pd_idx_ang= pd.DataFrame(idx_harm_ang, columns=['a1','a2','a3'])
               pd_ang_nat = pd.DataFrame(nat_ang, columns=['natang'])
               pd_coef_ang = pd.DataFrame(coef_ang, columns=['k'])
               self.protein_harmonic_angles = pd.concat([pd_idx_ang,pd_ang_nat,pd_coef_ang],axis=1)

     def aicg13_ang_array_to_pd(self):
          """
          To make aicg13 angles array to tabular format

          """
          if hasattr(self,'protein_aicg13_angles'):
               idx_ang = self.protein_aicg13_angles[0][:,6:9]
               epsilon = self.protein_aicg13_angles[1][:,3]*_kcal_to_kj
               r0 = self.protein_aicg13_angles[1][:,0] * _A_to_nm
               width = self.protein_aicg13_angles[1][:,4]* _A_to_nm
               pd_idx_ang = pd.DataFrame(idx_ang,columns=['a1', 'a2', 'a3'])
               pd_epsilon = pd.DataFrame(epsilon,columns=['epsilon'])
               pd_r0 = pd.DataFrame(r0, columns=['r0'])
               pd_width = pd.DataFrame(width, columns=['width'])
               self.protein_aicg13_angles = pd.concat([pd_idx_ang,pd_epsilon,pd_r0,pd_width],axis=1)

     def native_dihd_array_to_pd(self):
          """
          To make native dihedral angles array to tabular format

          """       
          if hasattr(self,'protein_native_dihd'):   
               idx_dihd = self.protein_native_dihd[0][:,7:11]
               nat_dihd = self.protein_native_dihd[1][:, 0]*np.pi/180
               coef_dihd = self.protein_native_dihd[1][:, 3:5]*_kcal_to_kj
               pd_idx_dihd = pd.DataFrame(idx_dihd,columns=['a1','a2','a3','a4'])
               pd_nat_dihd = pd.DataFrame(nat_dihd,columns=['natdihd'])
               pd_coef_dihd = pd.DataFrame(coef_dihd,columns=['k_dihd1','k_dihd3'])
               self.protein_native_dihd = pd.concat([pd_idx_dihd,pd_nat_dihd,pd_coef_dihd],axis=1)
     
     def aicg_dihd_array_to_pd(self):
          """
          To make aicg dihedral angles array to tabular format

          """
          if hasattr(self,'protein_aicg_dihd'):
               idx_dihd = self.protein_aicg_dihd[0][:,7:11]
               nat_dihd = self.protein_aicg_dihd[1][:,0]*np.pi/180
               epsilon = self.protein_aicg_dihd[1][:,3]*_kcal_to_kj
               width = self.protein_aicg_dihd[1][:,4]
               pd_idx_dihd = pd.DataFrame(idx_dihd,columns=['a1','a2','a3','a4'])
               pd_epsilon = pd.DataFrame(epsilon,columns=['epsilon'])
               pd_nat_dihd = pd.DataFrame(nat_dihd, columns=['natdihd'])
               pd_width = pd.DataFrame(width, columns=['width'])
               self.protein_aicg_dihd = pd.concat([pd_idx_dihd, pd_epsilon, pd_nat_dihd, pd_width], axis=1)
     
     def native_pair_array_to_pd(self):
          """
          To make native pairs array to tabular format

          """
          if hasattr(self,'protein_native_pair'):
               cont_para = self.protein_native_pair
               # index 
               idx_intra_cont_idx = np.argwhere(cont_para[0][:,2]==0)
               idx_inter_cont_idx = np.argwhere(cont_para[0][:,2]==1)
               idx_intra_cont = cont_para[0][idx_intra_cont_idx[:,0],:][:,5:7]
               idx_inter_cont = cont_para[0][idx_inter_cont_idx[:,0],:][:,5:7]

               # para
               epsilon_intra_cont = cont_para[1][idx_intra_cont_idx[:,0],3] * _kcal_to_kj
               r0_intra_con = cont_para[1][idx_intra_cont_idx[:,0],0] * _A_to_nm
               epsilon_inter_cont = cont_para[1][idx_inter_cont_idx[:,0],3] * _kcal_to_kj
               r0_inter_con = cont_para[1][idx_inter_cont_idx[:,0],0] * _A_to_nm
               para_intra_cont = np.stack((epsilon_intra_cont,r0_intra_con),axis=-1)
               para_inter_cont = np.stack((epsilon_inter_cont,r0_inter_con),axis=-1)

               # array to pandas
               pd_idx_intra_cont = pd.DataFrame(idx_intra_cont,columns=['a1','a2'])
               pd_para_intra_cont = pd.DataFrame(para_intra_cont,columns=['epsilon','sigma'])
               pd_idx_inter_cont = pd.DataFrame(idx_inter_cont,columns=['a1','a2'])
               pd_para_inter_cont = pd.DataFrame(para_inter_cont,columns=['epsilon','sigma'])
               
               # reset protein native pair
               if len(idx_intra_cont) != 0:
                  self.protein_intra_contact = pd.concat([pd_idx_intra_cont, pd_para_intra_cont],axis=1)
               if len(idx_inter_cont) != 0:
                  self.protein_inter_contact = pd.concat([pd_idx_inter_cont, pd_para_inter_cont],axis=1)
     
     def get_ninfo(self,ninfo_file_path):
          """
          get the native information in tabular format
          """
          self.parser_ninfo(ninfo_file_path)
          self.bonds_array_to_pd()
          self.harm_ang_array_to_pd()
          self.aicg13_ang_array_to_pd()
          self.native_dihd_array_to_pd()
          self.aicg_dihd_array_to_pd()
          self.native_pair_array_to_pd()
     
