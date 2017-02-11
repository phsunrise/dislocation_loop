'''
This script generates the positions of atoms that we wish to calculate
It is a list of (x,y,z,f) values for tier 1, or (x,y,z) for tier 2 and
above, where f is the factor of the atom that needs to be considered
'''
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from settings import NFILES, MAXTIER, preproc_dir, plot_dir
sample = 'W' 
from W_parameters import Vc, rot, crystaltype, a0, a2, ez

do_save = True 

R_cutoff = [300., 400.] # corresponding to cutoff for tier 1,2,...
maxind = int(np.ceil(R_cutoff[MAXTIER-1]*1./a0)) + 3 
                    # maximum index in any direction; 
                    # 3 is padding for possible shifting

if not crystaltype == 'BCC':
    print "Error: currently only works for BCC crystals"
    print "Aborting..."
    sys.exit(1)

if crystaltype == 'BCC':
    unit_cell = [np.array([0.,0.,0.]), \
                 np.array([0.5, 0.5, 0.5])]
    unit_cell_atoms = 2
elif crystaltype == 'FCC':
    unit_cell = [np.array([0.,0.,0.]), \
                 np.array([0.5, 0.5, 0.]), \
                 np.array([0.5, 0., 0.5]), \
                 np.array([0., 0.5, 0.5])]
    unit_cell_atoms = 4 

## first check if the preproc folder exists
if not os.path.isdir(preproc_dir):
    print "Creating folder %s ..."%preproc_dir
    os.system("mkdir %s"%preproc_dir)

for looptype in ['vac', 'int']:
    print "starting lootype %s" % looptype
    ## first calculate size of each file
    totallen1 = (4.*np.pi/3.*R_cutoff[0]**3)/Vc # this is only approximate
    sublistlen1 = int(np.ceil(totallen1*1./NFILES))
    xyz_list1 = []
    counter1 = 0
    i_file1 = 0
    totallen2 = (4.*np.pi/3.*(R_cutoff[1]**3-R_cutoff[0]**3))/Vc # this is only approximate
    sublistlen2 = int(np.ceil(totallen2*1./NFILES))
    xyz_list2 = []
    counter2 = 0
    i_file2 = 0
    for i_x in np.arange(-maxind, maxind+1):
        for i_y in np.arange(-maxind, maxind+1):
            for i_z in np.arange(-maxind, maxind+1):
                for i_atom in xrange(unit_cell_atoms): 
                    if looptype == 'vac':
                        xyz = a0*(np.array([i_x, i_y, i_z])+unit_cell[i_atom]) 
                    elif looptype == 'int':
                        # TODO: currently only applies to BCC with (110) plane
                        xyz = a0*(np.array([i_x-1./4, i_y+1./4, i_z])+unit_cell[i_atom]) 

                    xyz_p = rot.dot(xyz)
                    if np.linalg.norm(xyz_p) <= R_cutoff[0]:
                        xyz_list1.append(xyz_p)
                        counter1 += 1
                        if do_save:
                            if counter1 == sublistlen1:
                                np.save(preproc_dir+"%s_atoms_s_%s_pre_T1_%04d.npy"%(\
                                    sample, looptype, i_file1), xyz_list1)
                                print "saved T1 file %d, list length = %d" % (\
                                        i_file1, len(xyz_list1))
                                xyz_list1 = []
                                counter1 = 0
                                i_file1 += 1

                    elif R_cutoff[0] < np.linalg.norm(xyz_p) <= R_cutoff[1]:
                        xyz_list2.append(xyz_p)
                        counter2 += 1
                        if do_save:
                            if counter2 == sublistlen2:
                                np.save(preproc_dir+"%s_atoms_s_%s_pre_T2_%04d.npy"%(\
                                    sample, looptype, i_file2), xyz_list2)
                                print "saved T2 file %d, list length = %d" % (\
                                        i_file2, len(xyz_list2))
                                xyz_list2 = []
                                counter2 = 0
                                i_file2 += 1
    if do_save:
        if len(xyz_list1) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T1_%04d.npy"%(\
                sample, looptype, i_file1), xyz_list1)
            print "saved T1 file %d, list length = %d" % (\
                    i_file1, len(xyz_list1))
            i_file1 += 1
        print "done tier 1, total %d files" % (i_file1)
        if len(xyz_list2) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T2_%04d.npy"%(\
                sample, looptype, i_file2), xyz_list2)
            print "saved T2 file %d, list length = %d" % (\
                    i_file1, len(xyz_list2))
            i_file2 += 1
        print "done tier 2, total %d files" % (i_file2)

    ## plotting 
    if not do_save:
        print "length: %d, %d" % (\
                len(xyz_list1), len(xyz_list2))
        ### now expand the list to include negative values
        #def addfliplist(xyz_list):
        #    xyz_list = np.array(xyz_list)
        #    xyz_list_copy = np.copy(xyz_list)
        #    xyz_list_copy[:, 0:3] = -xyz_list_copy[:, 0:3]
        #    return np.vstack((xyz_list, xyz_list_copy))
        #xyz_list1 = addfliplist(xyz_list1)
        #xyz_list2 = addfliplist(xyz_list2)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for x,y,z in xyz_list1:
            if np.round(z/a2*2) in [1., 2.]:
                color = 'b'
                #ax.scatter(x/a2, y/a0*2, color=color, alpha=0.5)
                ax.scatter(x, y, color=color, alpha=0.5)
            elif np.round(z/a2*2) in [3., 4.]:
                color = 'r'
                #ax.scatter(x/a2, y/a0*2, color=color, alpha=0.5)
                ax.scatter(x, y, color=color, alpha=0.5)
        for x,y,z in xyz_list2:
            if np.round(z/a2*2) in [3., 4.]:
                #ax.scatter(x/a2, y/a0*2, color='g', alpha=0.5)
                ax.scatter(x, y, color='g', alpha=0.5)

        circ1 = plt.Circle((0, 0), R_cutoff[0], color='k', fill=False)
        circ2 = plt.Circle((0, 0), R_cutoff[1], color='k', fill=False)
        ax.add_artist(circ1)
        ax.add_artist(circ2)


        fig = plt.figure()
        ax = fig.add_subplot(111)
        for x,y,z in xyz_list1:
            if np.round(y/a0*2) == 0.:
                #ax.scatter(x/a2, z/a2, color='b')
                ax.scatter(x, z, color='b')
        for x,y,z in xyz_list2:
            if np.round(y/a0*2) == 0.:
                #ax.scatter(x/a2, z/a2, color='g', alpha=0.5)
                ax.scatter(x, z, color='g', alpha=0.5)

        circ1 = plt.Circle((0, 0), R_cutoff[0], color='k', fill=False)
        circ2 = plt.Circle((0, 0), R_cutoff[1], color='k', fill=False)
        ax.add_artist(circ1)
        ax.add_artist(circ2)

        plt.show()
