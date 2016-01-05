#!/usr/bin/env python
import numpy as np


######################################################
#
#	USER PARAMETERS
#
######################################################

#Dimensions of the Voxel datacube. The higher the number the higher the resolution.
#You can pick almost any resolution but don't go crazy. Take into account that the datacube will look empty if you have too many empty voxels in the regions where the particles are accumulated.
#Blender cannot read an infinite ammounte of voxels, so take the memory also in consideration.

nx=150
ny=150
nz=150

#Input and Output Files
Input="SPH.dat"
Output="SPH.bvox"


#If you want to set the data to logarithmic scale set logvalues=True	(Be sure your data does not have negative values)
logvalues=False


#If you know the approximate boundaries values of the data. (Recommended)
min_x,min_y,min_z=0,0,0
max_x,max_y,max_z=5,5,5

#Or if you want the code to get them for you set AutoMinMax=True. (Can produce be tricky to scale the blender cube later).
AutoMinMax=False


######################################################
#
#	SPH VOXELIZER
#
######################################################


#Load the columns
x,y,z,data= np.loadtxt(Input, unpack=True, usecols=(0,1,2,3))

#Calculate
if AutoMinMax:
	min_x,min_y,min_z=np.min(x),np.min(y),np.min(z)
	max_x,max_y,max_z=np.max(x),np.max(y),np.max(z)


#Scale (x,y,z) coordinates to be between 0 - (nx,ny,nz)
x= x - min_x
y= y - min_y
z= z - min_z
x= x * nx/(max_x-min_x)
y= y * ny/(max_y-min_y)
z= z * nz/(max_z-min_z)

x,y,z= x.astype(int), y.astype(int), z.astype(int)
#Normalize the data to be between 0-1

if logvalues:
	data=np.log10(data)


num=data.size 				#Get the lenght of the data array
vdata=np.zeros(nx*ny*nz)		#Array where we will store the voxel data
point_count= np.zeros(nx*ny*nz)		#Array where we count how main particles are inside each voxel. (For normalization purposes).


#Here we fill the vdata array and count the particles
for i in range(num):
	xx,yy,zz= x[i], y[i], z[i]
	vdata[xx+yy*nx+zz*nx*ny]+= data[i]
	point_count[xx+yy*nx+zz*nx*ny]+=1


point_count[np.where(point_count<1)[0]]=1	#To avoid divisions by zero
vdata=np.divide(vdata,point_count)		#Normalization to keep vdata values between 0-1


vdata= vdata - np.min(vdata)
vdata=vdata/np.max(vdata)


#Header of the blender voxel file. This is how blender knows what are the dimensions of the data cube.
header = np.array([nx,ny,nz,1])

#open and write to file
binfile = open(Output,'wb')
header.astype('<i4').tofile(binfile)
vdata.astype('<f4').tofile(binfile)
binfile.close()
