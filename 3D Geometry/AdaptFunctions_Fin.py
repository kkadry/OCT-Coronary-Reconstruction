# -*- coding: utf-8 -*-


# =============================================================================
# This module consists of a set of functions that will be used to
# automatically generate the required 3D geometry and computational mesh
# of atherosclerotic plaque.
# the main module used is gmsh which has functions to create a geometry and 
# correspondingly mesh it.
# it consists of the following functions:
#     
# 
#
#
# importdata: This function imports the mesh data of an abaqus odb file 
# (node positions+ indices, element indices and tags) along with
# the centroidal stress for each element in the form of txt files.
# The data is then used to inform the stress dependent remeshing algorithm.
#
# importdataV2: Same as importdata (in addition to strain data)
# but draws the same information from an npz file.
#
# Addcentroidpoints: This function takes as input the imported data (eltagsori,
# nodexyz,elemstress) from the abaqus odb file and creates reference points 
# at each element centroid with a stress higher than a specified value (filterstress)
#
# AddcentroidpointsV2: Similar to Addcentroidpoints but now filters points between
# a maximum and minimum stress (filterstressmin, filterstressmax) and Z position
# value (maxZ, minZ) 
#
#
# =============================================================================
import numpy as np
import gmsh
import sys
from scipy.spatial import cKDTree


def importdata(name):
    location=""
    nodeset=np.genfromtxt(location+"data_"+name+"_nodeset.txt",dtype="int")
    nodexyz=np.genfromtxt(location+"data_"+name+"_nodexyz.txt",delimiter=' ')
    elset=np.genfromtxt(location+"data_"+name+"_elset.txt",dtype="int")
    eltagsori=np.genfromtxt(location+"data_"+name+"_eltags.txt",dtype="int",delimiter=' ')
    elemstress=np.genfromtxt(location+"data_"+name+"_elemstress.txt")
    return(nodeset,nodexyz,elset,eltagsori,elemstress)

def importdata_V2(name):
    #Loading the matrices
    location=""
    data = np.load(location+'data_'+name+'.npz')
    nodeset=data["nodeset"]
    nodexyz=data["nodexyz"]
    elset=data["elset"]
    eltagsori=data["eltags"]
    elemstress=data["elemstress"]
    elemstrain=data["elemstrain"]
    return(nodeset,nodexyz,elset,eltagsori,elemstress,elemstrain)
    



def Addcentroidpoints(eltagsori,nodexyz,elemstress,filterstress):    
    #Filtering the low stress elements
    dex=np.where(elemstress>filterstress)[0] 
    eltags=eltagsori[dex]
    centroid=np.zeros((len(eltags),3))
    #Getting the centroids
    for i in range(len(eltags)):
        coordel=nodexyz[eltags[i]-1,:]
        centroid[i,:]=np.mean(coordel,axis=0)
    #Putting the points
    Npoint=np.shape(centroid)[0]
    #Declaring the characteristic mesh length
    #Declaring the point index list matching the points as given and their tags
    pointdex=np.empty(Npoint,dtype="int")
    
    #Adding all the points in coords (the -1 is to not put the duplicated point)
    for i in range(Npoint):
        pointdex[i]=gmsh.model.occ.addPoint(centroid[i,0], centroid[i,1], centroid[i,2])
    return(pointdex)
    
def AddcentroidpointsV2(eltagsori,nodexyz,elemstress,filterstressmin,filterstressmax,maxZ,minZ):    
    #Filtering the low stress elements
    dexmin=np.where(elemstress>filterstressmin)[0] 
    dexmax=np.where(elemstress<filterstressmax)[0]
    dex=np.intersect1d(dexmin,dexmax)
    eltags=eltagsori[dex]
    centroid=np.zeros((len(eltags),3))
    #Getting the centroids
    for i in range(len(eltags)):
        coordel=nodexyz[eltags[i]-1,:]
        centroid[i,:]=np.mean(coordel,axis=0)

    #Getting the centroids that are not at the caps
    centroidexmax=np.where(centroid[:,2]<maxZ)[0]
    centroidexmin=np.where(centroid[:,2]>minZ)[0]
    dex2=np.intersect1d(centroidexmin,centroidexmax)
    centroid=centroid[dex2]
    #Putting the points
    Npoint=np.shape(centroid)[0]
    #Declaring the characteristic mesh length
    #Declaring the point index list matching the points as given and their tags
    pointdex=np.empty(Npoint,dtype="int")
    
    #Adding all the points in coords (the -1 is to not put the duplicated point)
    for i in range(Npoint):
        pointdex[i]=gmsh.model.occ.addPoint(centroid[i,0], centroid[i,1], centroid[i,2])
    return(pointdex)

    


    
    



