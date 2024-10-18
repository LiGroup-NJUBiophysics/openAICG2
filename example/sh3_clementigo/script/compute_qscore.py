import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md


_kcal_to_kj = 4.1840
skipsnap = 500
T = [320,330,340,350,360,370]
native_contact = np.load('../output/native_contact.npy')
nat_con_idx = native_contact[:,:2]
sigma = native_contact[:,2]
metafile = open('metafile','w+')
for Ti in T:
    traj = md.load('../output/sh3_clementigo_%d.dcd'%Ti,top='../output/sh3_clementigo_%d.pdb'%Ti)
    distances = md.compute_distances(traj,atom_pairs=native_contact[:,:2])
    numsnapshot,_ = np.shape(distances)
    qscore = np.zeros(numsnapshot,dtype=float)
    for i in range(numsnapshot):
        qscore[i] = np.mean(1/(1+np.exp(50*(distances[i,:]-1.2*sigma))))
    #mdstep = np.arange(1,numsnapshot+1) 
    log = np.loadtxt('../output/sh3_clementigo_%d.log'%Ti,usecols=(1,2),dtype=float)
    mdstep = log[skipsnap:,0]
    potential = log[skipsnap:,1] / _kcal_to_kj
    qscore = qscore[skipsnap:]
    data = np.stack((mdstep[:],qscore[:],potential[:]),axis=-1)
    np.savetxt('../output/sh3_clementigo_%d.txt'%Ti,data)
    fig,ax = plt.subplots(figsize=(4.6,4))
    ax.plot(mdstep,qscore,'-k')
    ax.set_xlabel('MD step $10^3$',fontsize=24)
    ax.set_ylabel('Q', fontsize=24)
    plt.savefig('../output/Q_vs_mdstep_%d.png'%Ti,dpi=300,bbox_inches='tight')
    metafile.write('../output/sh3_clementigo_%d.txt %d')
    plt.close()

