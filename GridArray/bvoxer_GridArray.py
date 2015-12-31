#!/usr/bin/env python
import numpy as np


######################################################
#
#	USER PARAMETERS
#
######################################################

#Dimensions of the Voxel datacube. For this kind of simulations they should match the dimensions of the original datacube.
#Scientific data usually have a parameter file that (or include a header in the data file) which can help you to figure out these values.

nx=400
ny=400
nz=50

#Input and Output Files
Input="GridArray.dat"
Output="GridArray.bvox"

#If you want to set the data to logarithmic scale set logvalues=True	(Be sure your data does not have negative values)
logvalues=False


######################################################
#
#	GRID VOXELIZER
#
######################################################


#Load the data
data= np.loadtxt(Input)


if logvalues:
	data=np.log10(data)

#Normalize the data to be between 0-1
data= data- np.min(data)
data=data/np.max(data)

#Header of the blender voxel file. This is how blender knows what are the dimensions of the data cube.
header = np.array([nx,ny,nz,1])

#open and write to file
binfile = open(Output,'wb')
header.astype('<i4').tofile(binfile)
data.astype('<f4').tofile(binfile)
binfile.close()

