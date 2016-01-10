import numpy as np

Pi=3.14

def CustomFunction(x,y,z):
	#Replace the following line with CustomData=YourFunction.
	CustomData=np.sin(x)+ np.sin(y) + np.sin(z)
	return CustomData


###########################
#
#   USER PARAMETERS
#
###########################


nx,ny,nz =100,100,100

min_x,min_y,min_z= -Pi,-Pi,-Pi
max_x,max_y,max_z= Pi,Pi,Pi

TotalTime=200
TimeDomain=2.0*Pi
OutputFolder="VoxelOutputs/"


###########################
#
#   Generate VoxelData 
#
###########################

header=np.array([nx,ny,nz,1])

x,y,z=[],[],[]
for k in range(nz):
    zz=min_z + (max_z-min_z)*float(k)/float(nz)
    for j in range(ny):
        yy=min_y + (max_y-min_y)*float(j)/float(ny)
        for i in range(nx):
            xx=min_x + (max_x-min_x)*float(i)/float(nx)
            z.append(zz)
            y.append(yy)
            x.append(xx)
x,y,z= np.array(x),np.array(y),np.array(z)


for t in range (TotalTime):
    dr= TimeDomain* float(t)/float(TotalTime)
    data=CustomFunction(x+dr,y+dr,z+dr)
    data= (data- np.min(data))/(np.max(data)-np.min(data))

    #Open binary file. Write Header. Write Data.
    binfile = open(OutputFolder+"Voxel"+str(t)+".bvox",'wb')
    header.astype('<i4').tofile(binfile)
    data.astype('<f4').tofile(binfile)
    binfile.close()



    
