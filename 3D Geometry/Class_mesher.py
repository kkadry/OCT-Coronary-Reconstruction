# -*- coding: utf-8 -*-

# =============================================================================
# This module consists of a set of functions that will be used to
# automatically generate the required 3D geometry and computational mesh
# of atherosclerotic plaque.
# the main module used is gmsh which has functions to create a geometry and 
# correspondingly mesh it.
# it consists of the following classes and functions:
#
# plaque: the main geometric class representing a single atherosclerotic lesion
#
# importlip: this function imports the lipid step file
# importcalc: this function imports the calciumstep file
# importfibrous: this function imports the fibrous step file (unused)
# importartery: this function imports the artery step file
# importgeom: this function imports calcium,lipid, and artery
# 
# calciumcut: this function conducts a boolean subtraction operation between all calcific components and
# all lipidic components 
# fibrouscut: this function conducts a boolean subtraction operation between all fibrotic components and
# all lipidic components 
#
# cutbox: this function conducts a boolean cut operation between the artery and a rectangular box
# in order to take a submodel. Then the function conducts a boolean intersection
# operation between all plaque components and the submodel artery
# cutfull: this function conducts a boolean intersection operation between all plaque components
# and the full artery artery
#
#
# idxmap: this function maps the new geometric ID's after boolean operations 
# to the old indexes for lipid, calcium, and artery
#
# assignphys: this function assigns physical groups to each component for export to abaqus
#
# basemeshfield: this function implements default meshing options for gmsh 
#
# importdata: #This function imports stress and strain data for  abackground mesh to 
# inform mesh size fields 
# 
# meshfieldstress: This function applies a field for all elements in the background mesh
# between 2 values of stress based on the previous stressfield
#
# layermesh: This function applies meshfieldstress for multiple values of stress
# where higher stress values cause higher levels of mesh refinement
#
# meshgen: generates the mesh
# 
# arteryfind: finds the geometrical solid index corresponding to the artery through
# finding the biggest bounding box of every solid
#
# grun: runs the fltk gui to display the model
#
# gwrite: writes the mesh and physical groups to an abaqus inp file for processing
#
# removefield: removes defined sizing fields 
#=============================================================================
import numpy as np
import gmsh
import sys
from datetime import datetime
from scipy.spatial import cKDTree
from AdaptFunctionsV2 import*
import os
import subprocess
import paramiko
from PushBullet import*

class plaque:
    #Dimensionality
    dim=3
    #Initializing the model with model attributes
    def __init__(self,mname,lipname,calcname,fibname,artname,globalname):
        #Model name
        self.mname=mname
        #Global model name
        self.globalname=globalname
        #Names of the lipid, calcium and artery groups
        self.lipname=lipname
        self.calcname=calcname
        self.artname=artname
        self.fibname=fibname
        
        
        #Marks to check if certain materials are imported
        self.calciumcheck=False
        self.lipidcheck=False
        self.arterycheck=False
        self.fibrouscheck=False
        #Initialize the Gmsh kernel
        gmsh.initialize()
        #Option for printing to terminal
        gmsh.option.setNumber("General.Terminal", 0)
        #Setting the parallel option for boolean operators on
        gmsh.option.setNumber("Geometry.OCCParallel", 1)
        #adds the model
        #Gets a list of all models
        mlist=gmsh.model.list()
        #if modelname already exists, delete it
        if self.mname in mlist:
            gmsh.model.setCurrent(mname)
            #delete current model
            gmsh.model.remove()
        gmsh.model.add(mname)
        print("Model "+mname+" created")
        
    #Functions to import Lipid seperately
    def importlip(self):        
        gmsh.model.setCurrent(self.mname)
        self.lip = gmsh.model.occ.importShapes(self.lipname)#Lipid
        print("Lipid imported")  
        self.lipidcheck=True
        self.lip_ori=list(list(zip(*self.lip))[1])
        
    #Function to import calcium seperately  
    def importcalc(self):      
        gmsh.model.setCurrent(self.mname)
        self.calcium=gmsh.model.occ.importShapes(self.calcname)#Calcium
        print("Calcium imported")
        self.calciumcheck=True
        self.calcium_ori=list(list(zip(*self.calcium))[1])
    #Function to import fibrous seperately 
    def importfib(self):
        gmsh.model.setCurrent(self.mname)
        self.fibrous=gmsh.model.occ.importShapes(self.fibname)
        print("Fibrous imported")
        self.fibrouscheck=True
        self.fibrous_ori=list(list(zip(*self.fibrous))[1])
        
    #Function to import artery seperately    
    def importart(self):
        gmsh.model.setCurrent(self.mname)
        self.lum=gmsh.model.occ.importShapes(self.artname)#Lumen
        print("Artery imported")
        self.arterycheck=True
        self.lum_ori=list(list(zip(*self.lum))[1])
    
    #Shortcut function in case i want to import all of thems
    def importgeom(self):
        #Merging the step files into the model
        gmsh.model.setCurrent(self.mname)
        self.importlip()
        self.importcalc()
        self.importart()
    
    #This function takes the calcium volumes and cuts them against the lipid volumes
    def calciumcut(self):
        print("Cutting the calcium against the lipid")
        gmsh.model.setCurrent(self.mname)
        self.cut_out, self.cut_out_map =gmsh.model.occ.cut(self.lip,self.calcium,removeTool=False)
    
    #This function takes the lipid volumes and cuts them against the fibrous volumes
    def fibrouscut(self):
        print("Cutting the lipid against the fibrous")
        gmsh.model.setCurrent(self.mname)
        self.cut_fib_out, self.cut_fib_out_map =gmsh.model.occ.cut(self.fibrous,self.lip,removeTool=False)
        print("Cutting the calcium against the fibrous")
        #self.cut_fib_out2, self.cut_fib_out_ma2p =gmsh.model.occ.cut(self.fibrous,self.calcium,removeTool=False)
        
        
    #This function takes in the geometry and takes a section of all volumes that fit within a bounding box
    def cutbox(self,Centerslice=18,Nslice=8,buffer=0.2):
        if hasattr(self, 'lum')==False:
            self.lum=[]
        print("Starting the cutbox function")
        gmsh.model.setCurrent(self.mname)
        #Declaring the Center and length values for the plaque
        Center=Centerslice*0.4 
        Lbox=Nslice*0.4*2 #Multiplying by 2 
        #Finding the the start value for the bbox
        StartL=Center-Lbox/2
        #Adding the cutting box
        box=gmsh.model.occ.addBox(-5,-5,StartL,10,10,Lbox,1000)
        
        #Setting the limits of refinement and storing the values for later
        self.maxZ=StartL+Lbox-buffer
        self.minZ=StartL+buffer

        #Cutting only the lumen to save time
        #If there is no artery, skip
        if len(self.lum)>0:
            print("Cutting the lumen with the box")
            self.box_cut,self.box_cut_out=gmsh.model.occ.intersect([(3,box)],self.lum,removeObject=False)
        #Getting all entities to cut 
        gmsh.model.occ.synchronize()
        vols=gmsh.model.occ.getEntities(3)
        vols.remove((3,box))
        #Checking again if lumen exists  \
        if len(self.lum)>0:
            print("Cutting the volumes with the lumen")
            self.lum_int,self.lum_int_out=gmsh.model.occ.intersect(self.lum,vols,removeObject=False)
        
        vols=gmsh.model.occ.getEntities(3)
        vols.remove((3,box))

        #Cutting everything with the box
        self.box_int,self.box_int_out=gmsh.model.occ.intersect([(3,box)],vols,removeObject=True)
        
        #Removing all duplicate surfaces
        gmsh.model.occ.removeAllDuplicates()
    #This function takes in the geometry and cuts all the volumes against the lumen
    def cutfull(self):
        gmsh.model.setCurrent(self.mname)
        print("Intersecting with Main Artery")
        plaque=gmsh.model.occ.getEntities(3)
        artgroup=self.arteryfind()
        plaque.remove((3,artgroup))
        self.lum_int,self.lum_int_out=gmsh.model.occ.intersect([(3,artgroup)],plaque,removeObject=False)
        #Removing all duplicate surfaces
        gmsh.model.occ.removeAllDuplicates()
            
    
    #Function to unwrap the mappings
    def idxmap(self,mat_map):
        matgroup=[]
        #Unrolling each list
        for i in range(len(mat_map)):   
            #Check if list is empty
            if len(mat_map[i])>0:
                #If list has elements, loop over all elements and take second tag
                for j in range(len(mat_map[i])):
                        matgroup.append(mat_map[i][j][1])
        return(matgroup)
    
    def assign_phys(self):
        gmsh.model.setCurrent(self.mname)
        dim=self.dim
        #Physical group assignment functions
        lum_map_array=np.array(self.lum_int_out)
        #Getting the newly updated lists
        lip_map=lum_map_array[self.lip_ori]
        if self.calciumcheck:calcium_map=lum_map_array[self.calcium_ori]
        artery_map=lum_map_array[self.lum_ori]
        
        #getting the physical groups (still overlapping)
        lipgroup=self.idxmap(lip_map)
        if self.calciumcheck:calgroup=self.idxmap(calcium_map)
        artgroup=self.idxmap(artery_map)
        
        #Removing the intersections
        #Calcium should take any intersecting volume with lipid
        if self.calciumcheck:lipgroup=np.setdiff1d(lipgroup,calgroup)
        #Calcium and lipid should take any intersecting volumes with artyer
        if self.calciumcheck:artgroup=np.setdiff1d(artgroup,calgroup)
        artgroup=np.setdiff1d(artgroup,lipgroup)
        #Assigning the physical groups
        
        m1=gmsh.model.addPhysicalGroup(dim,lipgroup,tag=1)
        gmsh.model.setPhysicalName(dim, m1, "Lipid")
        
        if self.calciumcheck:
            m2=gmsh.model.addPhysicalGroup(dim,calgroup,tag=2)
            gmsh.model.setPhysicalName(self.dim, m2, "Calcium")
        
        m3=gmsh.model.addPhysicalGroup(dim,artgroup,tag=3)
        gmsh.model.setPhysicalName(self.dim, m3, "Artery")
        
        #m4=gmsh.model.addPhysicalGroup(dim,fibgroup,tag=4)
        #gmsh.model.setPhysicalName(self.dim, m4, "Fibrous")

    #This function implements the basic level of meshing
    #Meshing options temp
    def basemeshfield(self,lc2=0.1,ref1=2,minfac=1,maxfac=1):
        gmsh.model.setCurrent(self.mname)
        self.lc2=lc2
        self.removefield()
        print("Mesh Field Implementation Started") 
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0)
        #Setting the algorithm to Mesh Adapt
        gmsh.option.setNumber("Mesh.Algorithm",1)
        
        gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1)

        
        gmsh.model.occ.synchronize()
        #Finding the boundaries and artery index
        plaque=gmsh.model.getEntities(3)
        artgroup=self.arteryfind()
        #Removing the artery index
        plaque.remove((3,artgroup))
        #Finding all the plaque surfaces
        stag=gmsh.model.getBoundary(plaque,combined=False)
        stags=list(zip(*stag))[1]

                
        
        #removing the artery walls
        #Meshing options
        #Add a field that determines mesh size for the plaque component surfaces
        gmsh.model.mesh.field.add("Distance", 20)
        gmsh.model.mesh.field.setNumber(20, "NNodesByEdge", 100)
        gmsh.model.mesh.field.setNumbers(20, "NodesList", [])
        gmsh.model.mesh.field.setNumbers(20, "EdgesList", [])
        gmsh.model.mesh.field.setNumbers(20, "FacesList",np.array(stags))
        
        
        #Thresholding the distance field to allow for controlled mesh size gradient
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "IField", 20)
        gmsh.model.mesh.field.setNumber(2, "LcMin", lc2/ref1)
        gmsh.model.mesh.field.setNumber(2, "LcMax", lc2)
        gmsh.model.mesh.field.setNumber(2, "DistMin", lc2/ref1/2*minfac)
        gmsh.model.mesh.field.setNumber(2, "DistMax", lc2/ref1*maxfac)
        #Setting a background mesh in case this is final iteration
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
    
    #This function imports stress and strain data for  abackground mesh to inform size fields
    def importdata(self,name,location):
        #Loading the matrices      
        self.data = np.load(location+'data_'+name+'.npz')
        self.nodeset=self.data["nodeset"]
        self.nodexyz=self.data["nodexyz"]
        self.elset=self.data["elset"]
        self.eltagsori=self.data["eltags"]
        self.elemstress=self.data["elemstress"]
        self.elemstrain=self.data["elemstrain"]
        return(self.nodeset,self.nodexyz,self.elset,self.eltagsori,self.elemstress,self.elemstrain)
        
    #This function applies a field for all elements in the background mesh between 2 values of stress
    def meshfieldstress(self,filterstressmin,filterstressmax,refinelevel):
        gmsh.model.setCurrent(self.mname)
        
        #creating points of elements with a stress between min and max filterstress
        self.pointdex=AddcentroidpointsV2(self.eltagsori,self.nodexyz,self.elemstress,filterstressmin,filterstressmax,self.maxZ,self.minZ)
        #Synchronizing to register the points
        gmsh.model.occ.synchronize()
        #Creating a distance field and applying it to all the points produced by addcentroidpoints
        distfieldtag=gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumber(distfieldtag, "NNodesByEdge", 50)
        gmsh.model.mesh.field.setNumbers(distfieldtag, "NodesList",self.pointdex)
        gmsh.model.mesh.field.setNumbers(distfieldtag, "EdgesList", [])
        gmsh.model.mesh.field.setNumbers(distfieldtag, "FacesList",[])
        
        #Adding a threshold field to apply adequate mesh sizing
        threshtag=gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshtag, "IField", distfieldtag)
        gmsh.model.mesh.field.setNumber(threshtag, "LcMin", self.lc2/refinelevel)
        gmsh.model.mesh.field.setNumber(threshtag, "LcMax", self.lc2)
        gmsh.model.mesh.field.setNumber(threshtag, "DistMin", self.lc2/refinelevel)#/refinelevel)
        gmsh.model.mesh.field.setNumber(threshtag, "DistMax", self.lc2/2)#Original value of *2
        gmsh.model.mesh.field.setNumber(threshtag, "Sigmoid", 0)
        return(distfieldtag,threshtag)
    
    #This function adds multiple layered size fields
    def layermesh(self,filterstressvec,refvec,X):
        gmsh.model.setCurrent(self.mname)
        #Declaring the empty threshold tag vector
        threshvec=np.zeros(len(filterstressvec)-1,dtype='int')
        #Looping over filterstress vector in order to define several fields
        for i in range(len(filterstressvec)-1):
            distfieldtag,threshtag=self.meshfieldstress(filterstressvec[i],filterstressvec[i+1],refvec[i]*X)
            threshvec[i]=threshtag
        
        #Setting a background field which combines the base field with the layered fields
        mintag=gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(mintag, "FieldsList", np.append([2],threshvec))
        #Setting the background mesh to be the minimum mesh field
        gmsh.model.mesh.field.setAsBackgroundMesh(mintag)
        
        
        
    #Function that generates a 3d mesh
    def meshgen(self):
        gmsh.model.setCurrent(self.mname)
        print("Meshing Started")
        gmsh.model.mesh.generate(3)

        
    #This function checks the bbox volumes and returns the artery volume index (by finding the largest volume)
    def arteryfind(self):
        gmsh.model.setCurrent(self.mname)
        #Getting all the volumes
        vol=gmsh.model.occ.getEntities(3)
        #Declaring the volume and index vector for each volume
        bboxvol=np.zeros(len(vol))
        bboxidx=np.zeros(len(vol),dtype="int")
        
        #Looping over all volumes to find the bounding box volumes
        for i in range(len(vol)):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(3,vol[i][1])
            bboxvol[i]=(xmax-xmin)*(ymax-ymin)*(zmax-zmin)
            bboxidx[i]=vol[i][1]
            #Extracting the largest bounding box volume index (artery)
            artgroup=bboxidx[np.argmax(bboxvol)]
        return artgroup
    
    
    #This function runs the fltk gui
    def grun(self):
        gmsh.model.setCurrent(self.mname)
        gmsh.model.occ.synchronize()
        gmsh.fltk.run()
        
    #This function writes the mesh as an input file
    def gwrite(self,location):
        print("Writing the base input file")
        gmsh.model.setCurrent(self.mname)
        self.temploc=""
        gmsh.write(location+self.mname+".inp")
        
    #This function clears the mesh and any fields
    def removefield(self):
        gmsh.model.setCurrent(self.mname)
        gmsh.model.mesh.clear()
        gmsh.model.mesh.field.remove(7)
        gmsh.model.mesh.field.remove(4)
        gmsh.model.mesh.field.remove(3)
        gmsh.model.mesh.field.remove(22)
        gmsh.model.mesh.field.remove(2)
        gmsh.model.mesh.field.remove(20)




    
    


    
