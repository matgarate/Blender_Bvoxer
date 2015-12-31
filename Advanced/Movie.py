import bpy
import numpy as np

#i=0

for i in range(0,1001):
    bpy.context.scene.frame_current=i
    bpy.data.textures["Texture_Fargo"].voxel_data.filepath="gasdens"+ str(i) +".bvox"


    bpy.ops.render.render()    
    filepath_movie="Movie"+str(i)+".png"
    bpy.data.images[0].save_render(filepath_movie)