import os
InputFolder="GasData/"
OutputFolder="VoxelData/"



for i in range(638):
	print "Snap: "+str(i)
	command="./bvoxer -in "+InputFolder+"gasdens"+str(i)+".dat -out "+OutputFolder+"gasdens"+str(i)+".bvox -ns 512 -nr 256 -nc 64 -vnx 500 -vny 500 -vnz 100 -logr -logv -tri"
	os.system(command)
