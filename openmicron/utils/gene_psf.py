
def gene_psf(psf_file_name,top):
   """
   generate protein structure file for vmd 

   Parameters
   ----------
   psf_file_name: str
      The file name of PSF.

   top: topology
      The topology of OpenMM System.
   """
    
   atoms = list(top.atoms())
   bonds = list(top.bonds())
   space = ' '
   with open('%s.psf'%psf_file_name,'w') as psf:
         psf.write('PSF\n')
         # ATOM
         #psf.write(space*6)
         psf.write('%8d !NATOM\n'%len(atoms))
         for i in range(len(atoms)):
            psf.write('%8d%7d%7s%s%4s%5s%15d%14d%8d\n'%
                  (atoms[i].index+1,1,atoms[i].residue.name,
                  space*2,atoms[i].name,atoms[i].name,
                  0,100,0))
         psf.write('\n')
         psf.write('%8d !NBOND: bonds\n'%len(bonds))
         for i in range(len(bonds)):
            psf.write('%8d%8d'%(bonds[i][0].index+1,bonds[i][1].index+1))
            if (i+1)%4==0:
               psf.write('\n')
            elif i==int(len(bonds)-1):
               psf.write('\n')

