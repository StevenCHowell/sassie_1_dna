# $Id: survey_distances.py,v 1.2 2014-09-29 19:00:23 schowell Exp $
import os,sys,locale,string,numpy as np,math
import sassie.sasmol.sasmol as sasmol
import collision

def calc_dist(m1,coor):

    dist_array = []
    natoms = len(coor)
    try:
        for i in xrange(natoms-1):
            ix = coor[i][0]
            iy = coor[i][1]
            iz = coor[i][2]
            for j in xrange(i+1,natoms):
                fx = coor[j][0]
                fy = coor[j][1]
                fz = coor[j][2]
                dx2 = (ix-fx)*(ix-fx)
                dy2 = (iy-fy)*(iy-fy)
                dz2 = (iz-fz)*(iz-fz)
                dist_array.append(math.sqrt(dx2+dy2+dz2))
    except:
        print 'i,j = ',i,j
        print 'na = ',m1.natoms()
        sys.exit()

    return dist_array

def calc_dist2(coor1, coor2):

    dist_array = []
    natoms1 = len(coor1)
    natoms2 = len(coor2)

    try:
        for i in xrange(natoms1):
            x1 = coor1[i][0]
            y1 = coor1[i][1]
            z1 = coor1[i][2]
            for j in xrange(natoms2):
                x2 = coor2[j][0]
                y2 = coor2[j][1]
                z2 = coor2[j][2]
                dx2 = (x1-x2) ** 2
                dy2 = (y1-y2) ** 2
                dz2 = (z1-z2) ** 2
                dist_array.append(np.sqrt(dx2+dy2+dz2))
    except:
        print 'i,j = ',i,j
        print 'na1 = ',natoms1
        print 'na2 = ',natoms2
        sys.exit()

    return dist_array

def f_calc_dist2(coor1, coor2):
    import time

    dist = np.zeros(coor1.shape[0] * coor2.shape[0])
    dist = collision.overlap_dist(coor1, coor2, dist)
    timestr = time.strftime("%y%m%d_%H%M%S_") # save file prefix
    #np.savetxt(timestr+'dist.out', dist,fmt='%f',
    #           header='creted by survey_distances.f_calc_dist2: ' + timestr)
    return dist

def get_distances(m1,basis):

    frame = 0
    error,mask = m1.get_subset_mask(basis)
    if(len(error)):
        print 'subset error = ',error
    error,coor = m1.get_coor_using_mask(frame,mask)
    if(len(error)):
        print 'coor error = ',error

    return calc_dist(m1,coor[0])

def get_distances2(m1, basis1, basis2):

    frame = 0
    error,mask1 = m1.get_subset_mask(basis1)
    if(len(error)):
        print 'subset error = ',error

    error,mask2 = m1.get_subset_mask(basis2)
    if(len(error)):
        print 'subset error = ',error

    error,coor1 = m1.get_coor_using_mask(frame,mask1)
    if(len(error)):
        print 'coor error = ',error

    error,coor2 = m1.get_coor_using_mask(frame,mask2)
    if(len(error)):
        print 'coor error = ',error


    return f_calc_dist2(coor1[0], coor2[0])

def plot_historam(dist_array,title, N_bins):

    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    print 'making histogram plot'

    #nhist,nbins = np.histogram(dist_array, bins=np.arange(50), density=True)

    # the histogram of the data
    n, bins, patches = plt.hist(dist_array, N_bins, normed=True, color='green', alpha=0.75)

    plt.xlabel('Distance (Angstroms)')
    plt.ylabel('Probability')
    plt.title(title)
    #plt.xticks(np.arange(0.5, 4, 0.5))
    #plt.xlim([0.5, 4])
    #plt.autoscale(enable=True, axis='y', tight=None)
    plt.grid(True)

    plt.show()


if __name__ == '__main__':

    # pdbfile = 'new_dsDNA.pdb'       # no proteins
    # pdbfile = '1zbb_tetra_half.pdb' # not minimized
    # pdbfile = 'c11_withTails.pdb'   # wrong file (bad structure)
    pdbfile = 'new_c11_tetramer.pdb'
    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdbfile)

    #print 'CA results:'
    #basis = 'name[i] == "CA"'
    #ca_distances = get_distances(m1,basis)
    #print 'min distance ',min(ca_distances)
    #print 'max distance ',max(ca_distances)
    #plot_historam(ca_distances,'CA Distances')

    print 'DNA-Protein results:'
    basis1 = 'name[i] == "CA"' #protein CA atoms
    basis2 = 'name[i] == "P"'  #dna P atoms
    dna_protein_distances = get_distances2(m1, basis1, basis2)
    print 'min distance ',min(dna_protein_distances)
    print 'max distance ',max(dna_protein_distances)
    print 'pause'
    plot_historam(dna_protein_distances,'DNA-Protein Distances', len(dna_protein_distances)/100)

    #print 'Backbone results:'
    #basis = 'name[i] == "C" or name[i] == "O" or name[i] == "P"' #dna
    # #basis = 'name[i] == "N" or name[i] == "CA" or name[i] == "C"' #protein
    #backbone_distances = get_distances(m1,basis)
    #print 'min distance ',min(backbone_distances)
    #print 'max distance ',max(backbone_distances)
    #plot_historam(backbone_distances,'Backbone Distances', len(backbone_distances)/100)

    #print 'Heavy results:'
    #basis = 'name[i][0] != "H"'
    #heavy_distances = get_distances(m1,basis)
    #print 'min distance ',min(heavy_distances)
    #print 'max distance ',max(heavy_distances)
    #plot_historam(heavy_distances,'Heavy Distances', len(heavy_distances)/100)

    #print 'All atom results:'
    #basis = 'name[i] != "None"'
    #all_distances = get_distances(m1,basis)
    #print 'min distance ',min(all_distances)
    #print 'max distance ',max(all_distances)
    #plot_historam(all_distances,'All Distances', len(all_distances)/100)

