import os

for i in range(638):
	print "Snap: "+str(i)
	command="./bvoxer -in gasdens"+str(i)+".dat -out gasdens"+str(i)+".bvox -ns 512 -nr 256 -nc 64 -vnx 500 -vny 500 -vnz 100 -logr -logv"
	os.system(command)
