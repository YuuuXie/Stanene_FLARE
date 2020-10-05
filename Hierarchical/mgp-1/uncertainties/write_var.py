import numpy as np

mgp_v = np.load('mgp_v.npy')
f = open('traj_uncertainties.dat', 'w')
for v in range(mgp_v.shape[0]):
    f.write('\n\n')
    for a in range(mgp_v.shape[1]):
        f.write('{} {} {}\n'.format(mgp_v[v][a][0], 
                                    mgp_v[v][a][1], 
                                    mgp_v[v][a][2]))

f.close()                                    

