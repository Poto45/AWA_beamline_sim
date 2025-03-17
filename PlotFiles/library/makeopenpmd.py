import os
import numpy as np
from scipy import constants
from scipy.stats import halfnorm, norm
from openpmd_api import (Access, Dataset, Mesh_Record_Component, Series,
                         Unit_Dimension)

SCALAR = Mesh_Record_Component.SCALAR

def make_openPMD(x, y, z, px, py, pz, filename, qtotal, species_name='myparticle'):
    # open file for writing
    f = Series(
        filename,
        Access.create
    )
    # all required openPMD attributes will be set to reasonable default values
    # (all ones, all zeros, empty strings,...)
    # but one can also set domain-specific values
    f.meshes_path = "fields"
    f.particles_path = "particles"
    
    # new iteration
    cur_it = f.iterations[0]

    # particles
    electrons = cur_it.particles[species_name]
    electrons.set_attribute("comment", "What is love? Beam is love. Charged, focused love.")

    # define a charge                     quantity is "Current * Time" in ISQ
    electrons["charge"].unit_dimension = {Unit_Dimension.I: 1, Unit_Dimension.T: 1,}
    num = x.size
    print ("number of macroparticle:",  num)
    electrons["charge"][SCALAR].unit_SI = -1.0  # data below is already in Coulomb
    dset = Dataset(np.dtype("float64"), extent=[num])
    electrons["charge"][SCALAR].reset_dataset(dset)
    electrons["charge"][SCALAR].make_constant(1.60217662e-19)

    # now all other attributes
    weight = np.abs(qtotal / num / 1.60217662e-19)
    electrons["weighting"][SCALAR].reset_dataset(dset)
    electrons["weighting"][SCALAR].make_constant(weight)
    
    # mass
    electrons["mass"].unit_dimension = {Unit_Dimension.M: 1}
    electrons["mass"][SCALAR].unit_SI = 1.0
    electrons["mass"][SCALAR].reset_dataset(dset)
    electrons["mass"][SCALAR].make_constant(9.1093837015e-31)
    

    d = Dataset(x.dtype, extent=x.shape)
    electrons["position"].unit_dimension = {Unit_Dimension.L: 1}
    electrons["positionOffset"].unit_dimension = {Unit_Dimension.L: 1}

    electrons["position"]["x"].reset_dataset(d)
    electrons["position"]["y"].reset_dataset(d)
    electrons["position"]["z"].reset_dataset(d)

    electrons["positionOffset"]["x"].make_constant(0.)
    electrons["positionOffset"]["y"].make_constant(0.)
    electrons["positionOffset"]["z"].make_constant(0.)
    
    electrons["position"]["x"].store_chunk(x)
    electrons["position"]["y"].store_chunk(y)
    electrons["position"]["z"].store_chunk(z)
    
    d = Dataset(x.dtype, extent=x.shape)
    electrons["momentum"].unit_dimension = {Unit_Dimension.M: 1, Unit_Dimension.T: -1, Unit_Dimension.L: 1,}

    electrons["momentum"]["x"].reset_dataset(d)
    electrons["momentum"]["y"].reset_dataset(d)
    electrons["momentum"]["z"].reset_dataset(d)
    electrons["momentum"]["x"].store_chunk(px)
    electrons["momentum"]["y"].store_chunk(py)
    electrons["momentum"]["z"].store_chunk(pz)

    f.flush()

    del f
    
def split_openPMD(x, y, z, px, py, pz, filename, binlist, qtotal, species_name):
    
    mysum = 0
    fileprefix = os.path.split(filename)[-1]
    fnamelist = []
    prefixlist = []
    myqtotal = 0
    for i in range(binlist.size - 1):
        mask = (z <= binlist[i]) & (z > binlist[i+1])
        
        q_mask = mask.sum() / x.size * qtotal
        make_openPMD(x[mask], y[mask], z[mask],
                     px[mask], py[mask], pz[mask],
                     filename + str(i).zfill(1) + '.h5',
                     q_mask, species_name=species_name + str(i).zfill(1))
        
        fnamelist.append(filename + str(i).zfill(1) + '.h5')
        prefixlist.append(fileprefix + str(i).zfill(1))
        
        mysum += mask.sum()
        myqtotal += q_mask
        
    print(x.size, mysum)
    print(qtotal, myqtotal)
    return (prefixlist, fnamelist)



