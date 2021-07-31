# -*- coding: utf-8 -*-
import gmsh
import numpy as np

# =============================================================================
# This module consists of a set of functions that will be used to
# automatically generate the required geometry of atherosclerotic plaque, the 
# main module used is gmsh which has functions to create a geometry and 
# correspondingly mesh it.
# it consists of the following functions:
#     
# CreateSplines: this function takes in an Nx3 numpy matrix of point coordinates
# along with a list of lists containing a list of point connectivities to create
# several splines according to the number of lists specified. It also returns
# a vector specifying the spline tags for each indexed spline. 
# (for example: ((1,2,3,1),(4,5,6,4))) this specifies that the function should
# create two closed splines, with the first spline passing through  
# 1,2,3, and 1.
#
# MatSplines: this function takes in a matlab array containing the splines
# for each specific material within a single frame and recursively applies
# the CreateSplines function to draw every component and stores the indexes
# in associated lists.
#
# CreateLoops: this function takes in a list of lists indicating how each 
# loop should be constructed in terms of it's constituting lines/splines,
# (for example: ((1,2),(3,4)) indicates that two loops should be created, with
# the first loop connecting splines 1 and 2 etc.)
# 
# CreateExtrude: this function takes in a list indicating each closed
# loop that should be extruded and a Z value indicating extrusion length
# (for example: (1,2),3 indicates that the two loops with tags 1 and 2 should be
# extruded for a length of 3 in the Z direction.
#
# CreateLoft: this function takes in a list of lists indicating which set of
# closed loops should be extruded into volumes. For example: ((1,2))
# indicates that curves 1 and 2 should be lofted together into one volume
#
# Createcut: this function takes in a list of volumes and implements a series
# of binary boolean cut operations in order to eliminate any intersections
# between parts, the order of the list is the order of cutting, thus specifying
# (1,2,3,4) means that volume one will eat volumes 2,3, and 4 if they are 
# intersection and that volume 2 do the same to volumes 3 and 4. If a volume
# is completely removed then the function detects that and adjust accordingly
# if a volume is cut in half then the subsequent volumes will not cut or be cut
# in any later boolean operations.

# CreateSurf: This function creates a set of surfaces, it takes in the loop
# indices of the surfaces and adds a planar surface for each loop index.
#
# Partupdate:this function loops over ID list and assigns the corresponding loopdex to the
# Master list (part3d), part3d1[4]=(2,3) means that the fourth part in material
# 1 is created by lofting together loops 2 and 3. It takes in 3 lists,
# loopdex contains all the loop indexes for the specific frame and 
# material. IDlist contains all the part ID's for each loop in a single 
# frame for each material. part3d is the master list which is appended each
# loop.
#
# Comptranslate:this function translates the subpart/component index into the equivalent 
# loop index, this function is useless when only one material is plotted
# but becomes extremely useful when multiple materials/lumens are plotted as
# the curve/spline indexes will not align with the loop indexes. Takes in
# the slicepart list for a certain material which is a list of loop indexes
# indexed by slice position. It also takes in an input list of subparts
# to translate. The element value of each subpart is the location index of
# slicepart_flat. Note, if only one list should be passed as subparts it
# needs to be in the form of subparts=((1,2,3),) so it can be indexed 
# correctly.
#
# PushLumen: this function takes in a list of vessel border points and enforces
# a minimum wall thickness. Each inner vessel point is matched to it's 
# corresponding outer vessel point and the distance between them is classified as
# the thickness. Then, the outer vessel point is pushed out to a miniumum distance
# away if the thickness is below a specified limit.
# =============================================================================
#%%Function for quick visualization
def grun():
    gmsh.model.occ.synchronize()
    gmsh.fltk.run()
#%%This version of the func is written for when all the coordinates are already
#In order and form a single closed loop
def CreateSpline(coords,lc,loopcheck,surfcheck):
    import copy
    
    #Finding the number of points in coords
    Npoint=np.shape(coords)[0]
    #Declaring the charachteristic mesh length
    #Declaring the point index list matching the points as given and their tags
    pointdex=np.zeros(Npoint-1,dtype="int")
    
    #Adding all the points in coords (the -1 is to not put the duplicated point)
    for i in range(Npoint-1):
        pointdex[i]=gmsh.model.occ.addPoint(coords[i,0], coords[i,1], coords[i,2],lc)
    
    #Code to check clockwise/anticlockwise
    #Check if the 2nd point is higher than the first, if true then flip
    if (coords[1,1]-coords[0,1])>0:
        pointdex=np.flip(pointdex)
    #Declaring the number of splines in splineconnect
    Nspline=1
    

    
    #Declaring the index list matching the subparts/splines as given and tags
    splinedex=np.zeros(Nspline,dtype="int")
    loopdex=np.zeros(Nspline,dtype="int")
    surfdex=np.zeros(Nspline,dtype="int")
    #Copying the pointdex vector for modification
    splineconnect=copy.deepcopy(pointdex)

    

    

    #Closing the loop
    splineconnect=np.append(splineconnect,[splineconnect[0]])
    
    #Adding the spline connecting all the coords
    splinedex[0]=gmsh.model.occ.addSpline(splineconnect)
    
    
    #Checking if a loop has to be created from a closed spline
    if loopcheck==True:
            #Adding a closed curveloop with one spline per loop
        loopdex[0]=gmsh.model.occ.addCurveLoop([splinedex[0]])
            
    if surfcheck==True:
            #Creating a surface (if needed)
        surfdex[0]=gmsh.model.occ.addPlaneSurface([loopdex[0]])
    return(pointdex,splinedex,loopdex,surfdex)    
#%%Function that plots the splines associated with the components related to Mat X
def MatSpline(Framemat,lc,surfcheck,compcount,f):
    #Declaring a list of surfaces/components and sizes to return
    complist=[]
    sizelist=[]
    partlist=[]
    looplist=[]
    IDlist=[]        
    
    #Declaring the sampling rate
    samplemat=10
    if f in []:samplemat=2
    #Checking if the number of components is exactly 1
    if type(Framemat)!=np.ndarray:
        #Getting the boundary points arranged in clockwise order
        Matpts=Framemat.vertices

        #Checking if Matpts is larger than a certain limit
        if len(Matpts)>50:
            #Subsampling the pts to reduce complexity, use the last point to close the spline
            ptcap=Matpts[-1]
            Matpts=Matpts[0::samplemat]
            Matpts[-1]=ptcap

            
        #If surfcheck is switched on, draw surfaces
        if surfcheck==1:
            #Creating the spline using vertices and char length
            pointdex,Mdex,ldex,Sdex=CreateSpline(Matpts,lc,loopcheck=True,surfcheck=True)
        else:
            pointdex,Mdex,ldex,Sdex=CreateSpline(Matpts,lc,loopcheck=True,surfcheck=False)
        
        #If extrudecheck is switched on, store the surface index, else 
        #store the curve index
        if surfcheck==1:
            #Storing the component index
            complist.append(Sdex[0])
        else:
            complist.append(Mdex[0])
        #Storing the loop index
        looplist.append(ldex[0])
        #Storing the component size index (number of points as surrogate)
        sizelist.append(Framemat.vertices.shape[0])
        #Storing the part ID of the subpart
        IDlist.append(Framemat.ID)
        
        #Check if a subpart was just added
        if IDlist:
            #Updating the current subpart index
            compcount=compcount+1
        
        
        
        
    #Checking if the number of components is above 1
    elif Framemat.size>1:
        #Looping over all the components in the array
        for i in range(Framemat.size):
            #Getting the boundary points arranged in clockwise order
            Matpts=Framemat[i].vertices
            
            if len(Matpts)>30:
                #Subsampling the pts to reduce complexity
                ptcap=Matpts[-1]
                Matpts=Matpts[0::samplemat]   
                Matpts[-1]=ptcap  
                
         #Checking if the loop must be extruded   
            if surfcheck==1:
                #Creating the spline using vertices and char length
                pdex,Mdex,ldex,Sdex=CreateSpline(Matpts,lc,loopcheck=True,surfcheck=True)
            else:
                pdex,Mdex,ldex,Sdex=CreateSpline(Matpts,lc,loopcheck=True,surfcheck=False)
            
            if surfcheck==1:
                #Storing the component index
                complist.append(Sdex[0])
            else:
                complist.append(Mdex[0])
            #Storing the loop index
            looplist.append(ldex[0])
            #Storing the component size index (number of points as surrogate)
            sizelist.append(Framemat[i].vertices.shape[0])
            #Storing the part ID of the subpart
            IDlist.append(Framemat[i].ID)
            #Check if a subpart was just added
            if IDlist:
                #Updating the current subpart index
                compcount=compcount+1
    return(complist,sizelist,looplist,IDlist,compcount)
#%%
def CreateLoops(loopconnect):
    #Declaring number of loops to create
    Nloop=len(loopconnect)
    
    #Declaring the index list matching the closed loops as given and tags
    loopdex=np.zeros(Nloop,dtype="int")
    
    #Creating the loops
    for l in range(Nloop):
        loopdex[l]=gmsh.model.occ.addCurveLoop(np.array(loopconnect)[l])
    return(loopdex)
#%%
def CreateExtrude(extrudelist,Z):
    #Declaring the number of loops to extrude
    Nextrude=len(extrudelist)
    
    #Declaring the index list matching the volumes as given and tags
    voldex=np.zeros(Nextrude,dtype="int")
    #Extruding the volumes
    for e in range(Nextrude):
        exdim=gmsh.model.occ.extrude([(2, extrudelist[e])],dx=0,dy=0,dz=Z)
        
    #Storing the volume tag in a voldex array
        exdim=np.array(exdim)
        voldex[e]=exdim[exdim[:,0]==3,1]
    return(voldex)
#%%
def CreateLoft(loftlist):
    #Declaring the number of volumes to loft
    Nloft=len(loftlist)
    #Declaring the index list matching the volumes as given and tags
    loftdex=np.zeros(Nloft)
    #Lofting the parts
    for l in range(Nloft):
        loftdim=gmsh.model.occ.addThruSections(np.array(loftlist)[l],maxDegree=(2),makeRuled=False)
        loftdex[l]=loftdim[0][1]
    return(loftdex.astype(int),Nloft)

#%%
def CreateSurf(surflist):
    #Finding number of surfaces
    Nsurf=len(surflist)
    #Creating a surf index to return
    sdex=np.zeros(Nsurf,dtype=int)
    #Looping over all surfaces to create them
    for s in range(Nsurf):
        sdex[s]=gmsh.model.occ.addPlaneSurface([surflist[s]])
    return(sdex)

#%%
def CreateCut(cutlist,dim):
    while len(cutlist)>1:
        #Declaring the number of cuts that must be made
        Ncut=len(cutlist)
        #zipping the cutlist indexes into dimTag format: (2) becomes (3,2) etc
        cutTag=list(zip([dim]*Ncut,cutlist))
        
        #Getting the highest dominance dimtag to act as a tool that cuts all other vols
        toolTag=cutTag.pop(0)
        for i in range(len(cutTag)):
            #Cutting all the leftover volumes using the tooltag while not removing the tool
            gmsh.model.occ.cut([cutTag[i]],[toolTag],removeTool=False)
            
        #Calculating number of volumes to check if a part was completely cut
        gmsh.model.occ.synchronize()
        vols=gmsh.model.getEntities(dim)
        #Checking how many volumes
        Nvol=len(vols)
        if Nvol<len(cutlist): 
            #then one volume was completely cut and it must be fixed
            #unzipping the vols dimTags to compare with cutlist
            vols=list(zip(*vols))[1]
            
            #getting the completely cut out volumes
            cutvols=list(set(cutlist)-set(vols))     
            
            #removing the completely cut volumes from cutlist 
            for p in range(len(cutvols)):
                cutlist.remove(cutvols[p])
        
        #popping the list such that the most dominant tool is the next list element
        cutlist.pop(0)
#%%
def PartUpdate(loopdex,IDlist,part3d):
    #Loop over the entire IDlist and assign the loop index accordingly
    for l in range(len(IDlist)):
        part3d[IDlist[l]].append(loopdex[l])
    return(part3d)
#%%
def CompTranslate(subparts,slicepart):
    #Filter out the empty indices
    slicepart_filter = list(filter(None, slicepart))
    #Flatten the slicepart list of lists into a single list
    slicepart_flat = [val for sublist in slicepart_filter for val in sublist]
    #Convert to array since it's easier to play with indexes
    slicepart_array=np.array((slicepart_flat))
    
    #Loop over all the lists within subpart and translates the spline index
    #to loop index
    Nlist=len(subparts)
    #Initializing the translated loop list of lists that will be returned
    #Keep in mind that subpart index 1 should correspond to 0 indexing
    transloop=[]
    for i in range(Nlist):
        #Note: the -1 is to fix the indexing
        subdex=np.array(subparts[i])-1
        transloop.append(slicepart_array[subdex].tolist())
    return(transloop)

#%%
#This function takes in outer vessel and luminal points and returns the outer vessel 
#points with an enforced minimum wall thickness 
def PushLumen(Outerpts,Innerptsori):
    from scipy.spatial import cKDTree
    #Creating the cKDtree for nearest neighor queries
    tree = cKDTree(Innerptsori)
    #Pushing out the outer lumen points if thickness below certain value
    #Setting the minimum thickness value
    minthick=0.5
    #Subsampling every other point
    sample=15
    Outerpts=Outerpts[0::sample]
    #flattening the outer lumen points to get closest coordinates
    Outerptsflat=np.array(Outerpts).flatten()
    #Synchronizing so you can check closest coordinates
    gmsh.model.occ.synchronize()

    #Reshaping the flat array for innerspline
    Nout=int(len(Outerptsflat)/3)

    

    
    #Querying the cKDTree to get the closest pts to the outer vessel points
    thicknesstree,dex=tree.query(Outerpts)
    #Getting the distances and normalizing
    disttree=np.array(Outerpts-Innerptsori[dex])
    disttreenorm=disttree/thicknesstree[:,None]
    #Alternative way of getting the norm vector
    #Finding where the thickness is smaller than the minimum limit
    pushdextree=np.where(thicknesstree<minthick)
    pushdisttree=np.ones(Nout)*0
    pushdisttree[pushdextree]=np.ones(len(pushdextree))*minthick-thicknesstree[pushdextree]
    
    #Getting the centroid and the radial vector
    Centroid=np.mean(Outerpts,axis=0)
    disttree=Outerpts-Centroid
    norm=np.linalg.norm(disttree,axis=1)
    
    pushdistvectree=np.zeros((Nout,3))
    #Multiplying the push distance by the push direction
    for i in range(Nout):
        pushdistvectree[i,:]=pushdisttree[i]*disttree[i,:]/norm[i]
    

    #Pushing out all the Outerpts by a minimum thickness
    Outerptspush=Outerpts+pushdistvectree
    return(Outerptspush)