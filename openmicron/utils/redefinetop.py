class  RedefineTopology(object):
    def __init__(self):
        
        pass   

    def get_mini_component_set(self,chains):
        """
        To get minimal component set

        Parameters
        ----------
        top: topology
            The topology of OpenMM System.

        Return
        ------
        num_mini_set: int
            The number of minimal component set. 
        """
        sequence = []
        for _,i_chian in enumerate(chains):
            i_seq = ''
            for ires in i_chian._residues:
                i_seq+=ires.name
            sequence.append(i_seq)
        self.unique_sequence = list(set(sequence))
        len_mini_component = len(self.unique_sequence)
        self.num_mini_component_set = int(len(chains)/len_mini_component)
        self.num_atom_per_set = 0
        add_chain = []
        for i,iseq in enumerate(sequence):
            if iseq in self.unique_sequence and iseq not in add_chain:
               self.num_atom_per_set +=len(list(chains[i].atoms())) 
               add_chain.append(iseq)
        #return self.num_mini_component_set,self.num_atom_per_set
    
    def redefine_bond(self,top,bonds):
        """
        To redefine the bonds in topology.

        Parameters
        ----------
        top: topology
           The topology of OpenMM System.
        
        bonds: bond
           bonds in system.
        
        Return
        ------
        top: openmm topology
        """
        self.top = top
        self.atoms = list(self.top.atoms())
        self.top._bonds = []
        chains = self.top._chains
        self.get_mini_component_set(chains)
        for ic in range(self.num_mini_component_set):
            for i_bond in bonds:    
                index1 = ic*self.num_atom_per_set + i_bond[0].index
                index2 = ic*self.num_atom_per_set + i_bond[1].index
                self.top.addBond(self.atoms[index1],self.atoms[index2])

