# -*- coding: utf-8 -*-

import scipy.io as sio
import numpy as np
import gmsh
import copy
from Gmsh_Functions_Submission import*
from datetime import datetime
#For timing purposes
startTime = datetime.now()

#Initialize the Gmsh kernel
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
#Setting the parallel option for boolean operators on
gmsh.option.setNumber("Geometry.OCCParallel", 1)
    
    
gmsh.model.add("Model")

#Setting the characteristic length
lc=0.4

#Load in the patient materix, keep pythonic structure
pname=""
mat = sio.loadmat(pname,squeeze_me=True,struct_as_record=False)

#%%Options to be set here
#Vector containing any frames that must be skipped
frameskip=[] 

#Declaring a component list (useless if there are no surfaces)
#and sizelist for each frame 
compframe=[]
sizeframe=[]


#This variable controls where you produce 2D surfaces or 1D contours
surfswitch=1
#This variable controls whether you produce the inner and outer lumen splines
lumenswitch=1
#This variable controls whether the model is boolean fragmented
fragmentswitch=0

#Switch to tell the program to cut the components against one another
cutswitch=0
#Control variable to switch materials (calcium, fibrofatty,lipid, and fibrous) on and off
matcontrol=[1,0,0,0]

#Starting index for the slice
startdex=1
#Number of slices you loop over
Nslice=108

#%%
#Initializing any lists
#Stores all the spline indices of the inner lumens
Innerdexlist=[]
#Stores all the spline indices of the outer lumens
Outerdexlist=[]
#Stores all the subpart loop indices of the component for each slice 
#for each Mat
slicepart1=[]
slicepart2=[]
slicepart3=[]
slicepart4=[]

#Stores all the subpart spline indices for each slice for each mat 
#(same as slicepart)
compart1=[]
compart2=[]
compart3=[]
compart4=[]

#Maximum number of parts 
Nframes=70
#Stores all the components for each 3D part for each Mat
part3d1=[[] for i in range(Nframes)]
part3d2=[[] for i in range(Nframes)]
part3d3=[[] for i in range(Nframes)]
part3d4=[[] for i in range(Nframes)]

#Initializing the counters that keeps track of all the subparts, corresponds
#to current subpart index
comp1count=1
comp2count=1
comp3count=1
comp4count=1


#%%
#Load in the data regarding a specific frame "f"
print("Importing Materials")
for i in range(Nslice):
    f=i+startdex
    #Checking if frame should be skipped 
    if f in frameskip:
        print(print("Frame "+str(f)+" Skipped"))
        continue    
    
    Frame=mat["Frame"][f]
    print("Frame "+str(f)+" Imported")
    
    

    
    
    #Go over each material and draw the components for each frame
    if matcontrol[0]==1:
        comp1,size1,loop1,IDlist1,comp1count=MatSpline(Frame.Mat1,lc,surfcheck=surfswitch,compcount=comp1count,f=f)
        compframe=compframe+comp1
        sizeframe=sizeframe+size1
        slicepart1.append(loop1)
        compart1.append(comp1)
        
        #Updating the part indexed master list
        part3d1=PartUpdate(loop1,IDlist1,part3d1)
        
    if matcontrol[1]==1:
        comp2,size2,loop2,IDlist2,comp2count=MatSpline(Frame.Mat2,lc,surfcheck=surfswitch,compcount=comp2count,f=f)
        compframe=compframe+comp2
        sizeframe=sizeframe+size2
        slicepart2.append(loop2)
        compart2.append(comp2)
        
        #Updating the part indexed master list
        part3d2=PartUpdate(loop2,IDlist2,part3d2)
        
    if matcontrol[2]==1:    
        comp3,size3,loop3,IDlist3,comp3count=MatSpline(Frame.Mat3,lc,surfcheck=surfswitch,compcount=comp3count,f=f)
        compframe=compframe+comp3
        sizeframe=sizeframe+size3
        slicepart3.append(loop3)
        compart3.append(comp3)
        
        #Updating the part indexed master list
        part3d3=PartUpdate(loop3,IDlist3,part3d3)
        
    if matcontrol[3]==1:
        comp4,size4,loop4,IDlist4,comp4count=MatSpline(Frame.Mat4,lc,surfcheck=surfswitch,compcount=comp4count,f=f)
        compframe=compframe+comp4
        sizeframe=sizeframe+size4
        slicepart4.append(loop4)
        compart4.append(comp4)
        
        #Updating the part indexed master list
        part3d4=PartUpdate(loop4,IDlist4,part3d4)
     

    
    #This does some setup for the 2D cutting (surface against surface)
    if surfswitch==1 and cutswitch==1:
        #Sorting the 2D components into cutting order
        #Turning the component list into array for sorting
        comparray=np.array(compframe)
        #Getting the sorting key based on size of each component (small first)
        cutorder=np.argsort(sizeframe)
        #Getting the cutting order list
        CutList=comparray[cutorder].tolist()
        
        #Cutting all the plaque components to ensure no intersection
        CreateCut(CutList,2)
    
    #Checking for lumenswitch to draw the vessel borders
    if lumenswitch==1:
        
        #Drawing the lumen borders
        #Getting the inner lumen points
        Innerptsori=Frame.Inner.vertices
        #Declaring the sampling size
        sample=15
        #Subsampling the inner lumen points
        Innerpts=Innerptsori[0::sample]
        
        #Whether or not the inner lumens have a 2d surface depends on surfswitch 
        if surfswitch==1:
            #Drawing the spline and loop
            pdex,innerdex,innerldex,innersdex=CreateSpline(Innerpts,lc*3,loopcheck=True,surfcheck=True)  
        else:
            pdex,innerdex,innerldex,innersdex=CreateSpline(Innerpts,lc*3,loopcheck=True,surfcheck=False) 

        #Storing the innerldex (loop index) for the slice in a larger structure (for lofting)
        Innerdexlist.append(innerldex[0])
        
        
        
        #Drawing the outer vessel borders
        #Getting the Outer vessel points
        Outerpts=Frame.Outer.vertices
        #Pushing the outer lumen points
        Outerptspush=PushLumen(Outerpts,Innerptsori)
            
        
        if surfswitch==1:
            #Drawing the spline and loop
            pdex,outerdex,outerldex,outersdex=CreateSpline(Outerptspush,lc*3,loopcheck=True,surfcheck=True)
        else:
            pdex,outerdex,outerldex,outersdex=CreateSpline(Outerptspush,lc*3,loopcheck=True,surfcheck=False)
        #Storing the outerdex for the slice in a larger structure
        Outerdexlist.append(outerldex[0])
    




#Drawing the vessel borders at the final frame  
if lumenswitch==1:
    #Reaching ahead to get the inner and outer lumen for the final slice 
    Frame=mat["Frame"][f+1]
    #Getting the inner lumen points
    Innerptsori=Frame.Inner.vertices
    #Subsampling the iner lumen points
    Innerpts=Innerptsori[0::sample]
    #Drawing the spline and loop
    pdex,innerdex,innerldex,sdex=CreateSpline(Innerpts,lc*3,loopcheck=True,surfcheck=False)
    #Storing the innerdex for the slice in a larger structure
    Innerdexlist.append(innerldex[0])
    
    
    #Getting the Outer lumen points
    Outerpts=Frame.Outer.vertices
    #Pushing the outer lumen outward
    Outerptspush=PushLumen(Outerpts,Innerptsori)
    #Drawing the spline and loop
    pdex,outerdex,outerldex,sdex=CreateSpline(Outerptspush,lc*3,loopcheck=True,surfcheck=False)
    #Storing the outerdex for the slice in a larger structure
    Outerdexlist.append(outerldex[0])


print(datetime.now() - startTime)


if fragmentswitch==1:
    print("Starting the fragmentation...")
    #This command removes any duplicated entities and coheres the geometry
    gmsh.model.occ.removeAllDuplicates()
else:
    print("Fragmentation Skipped")


#Geometry display options
gmsh.option.setNumber("Geometry.Points", 0)



#This command synchronizes the gmsh object with the program kernel
gmsh.model.occ.synchronize()
#Run the gui to display
gmsh.fltk.run()
#Closes the model
gmsh.finalize()

