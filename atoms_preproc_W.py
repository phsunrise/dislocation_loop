'''
This script generates the positions of atoms that we wish to calculate
It is a list of (x,y,z,f) values for tier 1, or (x,y,z) for tier 2 and
above, where f is the factor of the atom that needs to be considered
'''
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from info import NFILES, MAXTIER 
from settings import preproc_dir, plot_dir

do_save = True 
tiers = [1, 2, 3]

size = 2**7 # a multiple of 2

sample = 'W' 
looptype = 'int'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *

## first check if the preproc folder exists
if not os.path.isdir(preproc_dir):
    print "Creating folder %s ..."%preproc_dir
    os.system("mkdir %s"%preproc_dir)

if 1 in tiers:
    ## Tier 1: single atoms
    ## first calculate size of each file
    totallen1 = (size-1)*(size*2-1)**2/2 + (8*size-1)**2/2 # second term is the plane of half atoms
    sublistlen1 = int(np.ceil(totallen1*1./NFILES))
    xyz_list1 = []
    counter1 = 0
    i_file1 = 0
    for z_p in np.arange(1, size):
        for xpy in np.arange(-size+1, size):
            for xmy in np.arange(-size+1, size):
                '''
                if z_p is even, we have plane A, which has an atom at x=y=0
                if z_p is odd, we have plane B, which has an atom at x=y=1/2
                '''
                if (xpy+xmy+z_p) % 2 == 0:
                    x_p = (xpy+xmy)/2.
                    y_p = (xpy-xmy)/2.
                    if looptype == 'vac':
                        xyz = x_p*ex_p + y_p*ey_p + z_p*ez_p
                    elif looptype == 'int':
                        xyz = (x_p-0.25)*ex_p + (y_p-0.25)*ey_p + (z_p-0.5)*ez_p

                    f = 1. # factor of atom that needs to be considered
                    if z_p in [1, size-1]:
                        f = f/2.
                    if abs(xpy) == size-1:
                        f = f/2.
                    if abs(xmy) == size-1:
                        f = f/2.
                    #if z_p == 1:
                    #    f += 0.5 # add back the half atoms on bottom plane

                    xyz_list1.append(np.append(xyz, f))
                    counter1 += 1

                    if do_save:
                        if counter1 == sublistlen1:
                            np.save(preproc_dir+"%s_atoms_s_%s_pre_T1_%04d.npy"%(\
                                sample, looptype, i_file1), xyz_list1)
                            print "saved file %d, list length = %d" % (\
                                    i_file1, len(xyz_list1))
                            xyz_list1 = []
                            counter1 = 0
                            i_file1 += 1

    # the plane of half atoms at bottom, up to MAXTIER
    z_p = 1 
    f = 0.5 # factor of atom that needs to be considered
    for xpy in np.arange(-2**(MAXTIER-1)*size+1, 2**(MAXTIER-1)*size):
        for xmy in np.arange(-2**(MAXTIER-1)*size+1, 2**(MAXTIER-1)*size):
            '''
            since z_p is odd, we have plane B, which has an atom at x=y=1/2
            '''
            if (xpy+xmy+z_p) % 2 == 0:
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                if looptype == 'vac':
                    xyz = x_p*ex_p + y_p*ey_p + z_p*ez_p
                elif looptype == 'int':
                    xyz = (x_p-0.25)*ex_p + (y_p-0.25)*ey_p + (z_p-0.5)*ez_p

                xyz_list1.append(np.append(xyz, f))
                counter1 += 1

                if do_save:
                    if counter1 == sublistlen1:
                        np.save(preproc_dir+"%s_atoms_s_%s_pre_T1_%04d.npy"%(\
                            sample, looptype, i_file1), xyz_list1)
                        print "saved file %d, list length = %d" % (\
                                i_file1, len(xyz_list1))
                        xyz_list1 = []
                        counter1 = 0
                        i_file1 += 1

    if do_save:
        if len(xyz_list1) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T1_%04d.npy"%(\
                sample, looptype, i_file1), xyz_list1)
            print "saved file %d, list length = %d" % (\
                    i_file1, len(xyz_list1))
            i_file1 += 1
        print "done tier 1, total %d files" % (i_file1)

if 2 in tiers:
    ## Tier 2: single cells
    totallen2 = (size-1)*(size*2-1)**2 - (size/2-1)*(size-1)**2
    sublistlen2 = int(np.ceil(totallen2*1./NFILES))
    xyz_list2 = []
    counter2 = 0
    i_file2 = 0
    ## first expand along z
    for z_p in np.arange(2, size*2, 2):
        for xpy in np.arange(-2*size+2, 2*size, 2):
            for xmy in np.arange(-2*size+2, 2*size, 2):
                if z_p < size and abs(xpy) < size and abs(xmy) < size:
                    continue
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                if looptype == 'vac':
                    xyz = x_p*ex_p + y_p*ey_p + z_p*ez_p
                elif looptype == 'int':
                    xyz = (x_p-0.25)*ex_p + (y_p-0.25)*ey_p + (z_p-0.5)*ez_p
                xyz_list2.append(xyz)
                counter2 += 1

                if do_save:
                    if counter2 == sublistlen2:
                        np.save(preproc_dir+"%s_atoms_s_%s_pre_T2_%04d.npy"%(\
                            sample, looptype, i_file2), xyz_list2)
                        print "saved file %d, list length = %d" % (\
                                i_file2, len(xyz_list2))
                        xyz_list2 = []
                        counter2 = 0
                        i_file2 += 1

    if do_save:
        if len(xyz_list2) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T2_%04d.npy"%(\
                sample, looptype, i_file2), xyz_list2)
            print "saved file %d, list length = %d" % (\
                    i_file2, len(xyz_list2))
            i_file2 += 1
        print "done tier 2, total %d files" % (i_file2)

if 3 in tiers:
    ## Tier 3: 2^3 cells
    totallen3 = (size-1)*(size*2-1)**2 - (size/2-1)*(size-1)**2
    sublistlen3 = int(np.ceil(totallen3*1./NFILES))
    xyz_list3 = []
    counter3 = 0
    i_file3 = 0
    ## first expand along z
    for z_p in np.arange(4, size*4, 4):
        for xpy in np.arange(-4*size+4, 4*size, 4):
            for xmy in np.arange(-4*size+4, 4*size, 4):
                if z_p < 2*size and abs(xpy) < 2*size and abs(xmy) < 2*size:
                    continue
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                if looptype == 'vac':
                    xyz = x_p*ex_p + y_p*ey_p + z_p*ez_p
                elif looptype == 'int':
                    xyz = (x_p-0.25)*ex_p + (y_p-0.25)*ey_p + (z_p-0.5)*ez_p
                xyz_list3.append(xyz)
                counter3 += 1

                if do_save:
                    if counter3 == sublistlen3:
                        np.save(preproc_dir+"%s_atoms_s_%s_pre_T3_%04d.npy"%(\
                            sample, looptype, i_file3), xyz_list3)
                        print "saved file %d, list length = %d" % (\
                                i_file3, len(xyz_list3))
                        xyz_list3 = []
                        counter3 = 0
                        i_file3 += 1

    if do_save:
        if len(xyz_list3) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T3_%04d.npy"%(\
                sample, looptype, i_file3), xyz_list3)
            print "saved file %d, list length = %d" % (\
                    i_file3, len(xyz_list3))
            i_file3 += 1
        print "done tier 3, total %d files" % (i_file3)

if 4 in tiers:
    ## Tier 4: 4^3 cells
    totallen4 = (size-1)*(size*2-1)**2 - (size/2-1)*(size-1)**2
    sublistlen4 = int(np.ceil(totallen4*1./NFILES))
    xyz_list4 = []
    counter4 = 0
    i_file4 = 0
    ## first expand along z
    for z_p in np.arange(8, size*8, 8):
        for xpy in np.arange(-8*size+8, 8*size, 8):
            for xmy in np.arange(-8*size+8, 8*size, 8):
                if z_p < 4*size and abs(xpy) < 4*size and abs(xmy) < 4*size:
                    continue
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                if looptype == 'vac':
                    xyz = x_p*ex_p + y_p*ey_p + z_p*ez_p
                elif looptype == 'int':
                    xyz = (x_p-0.25)*ex_p + (y_p-0.25)*ey_p + (z_p-0.5)*ez_p
                xyz_list4.append(xyz)
                counter4 += 1

                if do_save:
                    if counter4 == sublistlen4:
                        np.save(preproc_dir+"%s_atoms_s_%s_pre_T4_%04d.npy"%(\
                            sample, looptype, i_file4), xyz_list4)
                        print "saved file %d, list length = %d" % (\
                                i_file4, len(xyz_list4))
                        xyz_list4 = []
                        counter4 = 0
                        i_file4 += 1

    if do_save:
        if len(xyz_list4) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T4_%04d.npy"%(\
                sample, looptype, i_file4), xyz_list4)
            print "saved file %d, list length = %d" % (\
                    i_file4, len(xyz_list4))
            i_file4 += 1
        print "done tier 4, total %d files" % (i_file4)

if 5 in tiers:
    ## Tier 5: 8^3 cells
    totallen5 = (size-1)*(size*2-1)**2 - (size/2-1)*(size-1)**2
    sublistlen5 = int(np.ceil(totallen5*1./NFILES))
    xyz_list5 = []
    counter5 = 0
    i_file5 = 0
    ## first expand along z
    for z_p in np.arange(16, size*16, 16):
        for xpy in np.arange(-16*size+16, 16*size, 16):
            for xmy in np.arange(-16*size+16, 16*size, 16):
                if z_p < 8*size and abs(xpy) < 8*size and abs(xmy) < 8*size:
                    continue
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                if looptype == 'vac':
                    xyz = x_p*ex_p + y_p*ey_p + z_p*ez_p
                elif looptype == 'int':
                    xyz = (x_p-0.25)*ex_p + (y_p-0.25)*ey_p + (z_p-0.5)*ez_p
                xyz_list5.append(xyz)
                counter5 += 1

                if do_save:
                    if counter5 == sublistlen5:
                        np.save(preproc_dir+"%s_atoms_s_%s_pre_T5_%04d.npy"%(\
                            sample, looptype, i_file5), xyz_list5)
                        print "saved file %d, list length = %d" % (\
                                i_file5, len(xyz_list5))
                        xyz_list5 = []
                        counter5 = 0
                        i_file5 += 1

    if do_save:
        if len(xyz_list5) > 0:
            np.save(preproc_dir+"%s_atoms_s_%s_pre_T5_%04d.npy"%(\
                sample, looptype, i_file5), xyz_list5)
            print "saved file %d, list length = %d" % (\
                    i_file5, len(xyz_list5))
            i_file5 += 1
        print "done tier 5, total %d files" % (i_file5)


## plotting 
if not do_save:
    print "length: %d, %d, %d" % (\
            len(xyz_list1), len(xyz_list2), len(xyz_list3))
    ## now expand the list to include negative values
    def addfliplist(xyz_list):
        xyz_list = np.array(xyz_list)
        xyz_list_copy = np.copy(xyz_list)
        xyz_list_copy[:, 0:3] = -xyz_list_copy[:, 0:3]
        return np.vstack((xyz_list, xyz_list_copy))
    xyz_list1 = addfliplist(xyz_list1)
    xyz_list2 = addfliplist(xyz_list2)
    xyz_list3 = addfliplist(xyz_list3)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for x,y,z,f in xyz_list1:
        if np.round(z/a2*2) in [1., 2.]:
            ax.scatter(x/a2, y/a0*2, color='b', alpha=0.5)
            ax.text(x/a2, y/a0*2, "%d"%(1./f))
        elif np.round(z/a2*2) in [3., 4.]:
            ax.scatter(x/a2, y/a0*2, color='r', alpha=0.5)
            ax.text(x/a2, y/a0*2, "%d"%(1./f))
    for x,y,z in xyz_list2:
        if np.round(z/a2*2) in [3., 4.]:
            ax.scatter(x/a2, y/a0*2, color='g', alpha=0.5)
    for x,y,z in xyz_list3:
        if np.round(z/a2*2) in [7., 8.]:
            ax.scatter(x/a2, y/a0*2, color='m', alpha=0.5)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for x,y,z,f in xyz_list1:
        if np.round(y/a0*2) == 0.:
            ax.scatter(x/a2, z/a2, color='b')
            ax.text(x/a2, z/a2, "%d"%(1./f))
    for x,y,z in xyz_list2:
        if np.round(y/a0*2) == 0.:
            ax.scatter(x/a2, z/a2, color='g', alpha=0.5)
    for x,y,z in xyz_list3:
        if np.round(y/a0*2) == 0.:
            ax.scatter(x/a2, z/a2, color='m', alpha=0.5)

    plt.show()
