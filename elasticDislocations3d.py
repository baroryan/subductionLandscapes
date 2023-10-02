#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 09:46:10 2022

@author: bar
"""
import okada4py
#import time
import numpy as np
import  xarray as xr
import matplotlib.pyplot as plt
import elasticDislocations3d
import ComputeWedge
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import tri
import pygmsh
import pandas as pd
import miscFunctions
#%%
def ConvertDataToArrayFloat64(data):
        if np.isscalar(data):
            dataToReturn=np.array([np.float64(data)])
        elif len(data)==1 and np.isscalar(data) is False:
            dataToReturn=np.array([np.float64(data)])
        else:
            dataToReturn=(np.float64(data))
        
        return dataToReturn
        
        
    

#%%
class grid:
    """ this class get array of x and y and perpeares the grid to work with okada92.
    It also makes sure everything is array of float 92
    because okada92 is very fincky about it
    x and y are vectors
    
        |------------------------------->x
        |      
        |
        |        \
        |         \
        |          \ fault map view
        |           \ 
        |            \
        |             \
        |
        y (map view)    
        
      can find an exmaple at /Users/bar/GoogleDrive/phd/research/Arthur/computeUplift/variesCouplingToTest/tryDifferentElasticSoultion.py  """
    
    
    def __init__(self,x,y,z=0):
        
        self.x=x
        self.y=y
        self.z=np.float64(z)

        
        self.ComputeGrid()
        
    def ComputeGrid(self):
        self.MakeSureAllDataIsFloat()
        self.GenerateGrid(self.x,self.y,self.z)
    
    def MakeSureAllDataIsFloat(self):
        self.x=ConvertDataToArrayFloat64(self.x)
        self.y=ConvertDataToArrayFloat64(self.y)
        self.z=ConvertDataToArrayFloat64(self.z)
        
            
    def GenerateGrid(self,x,y,z):
         
         self.xs,self.ys=np.meshgrid(x,y)
         self.gridShape=self.xs.shape
         self.xsFlat = self.xs.flatten()
         self.ysFlat = self.ys.flatten()
         self.zsFlat = np.ones_like(self.ysFlat)*self.z
         
    def PlotGrid(self,ax=None,**args):
        if ax is None:
            fig,ax=plt.subplots()
            
        ax.scatter(self.xsFlat,self.zsFlat,**args)
        
        
    def GenerateEmptyXRGrid(self):
        zeroValue=np.zeros(self.gridShape)
        dz=xr.DataArray(data=zeroValue,dims=["y","x"],coords=[self.y,self.x])
        dx=xr.DataArray(data=zeroValue,dims=["y","x"],coords=[self.y,self.x])
        dy=xr.DataArray(data=zeroValue,dims=["y","x"],coords=[self.y,self.x])
        
        arrays=[]
        for dd in [dx,dy,dz]:
            array=xr.DataArray(data=dd,dims=["y","x"],coords=[self.y,self.x])
            array.x.attrs['units']='km'
            array.y.attrs['units']='km'
            array.attrs['units']='m'
            arrays.append(array)
            
        
        return xr.Dataset(data_vars={'dz':arrays[2],'dy':arrays[1],'dx':arrays[0]})
        
        
        
         


                 
#%%
class griddx(grid):
    """ this class get array of x/y start/end and dx/dy and perpeares the grid to work with okada92. 
    It also makes sure everything is array of float 92
    because okada92 is very fincky about it"""
    def __init__(self,startX,endX,dx,startY,endY,dy,z=0):
        self.x=np.arange(startX,endX,dx)
        self.y=np.arange(startY,endY,dy)
        self.z=np.float64(z)
            
        self.ComputeGrid()
#%%
class gridxz(grid):
    """ this classes makes everything ready for 2d stresses on depth and distance from fault stresses:
        |------------------------------->x
        |      
        |
        |        \
        |         \
        |          \ fault
        |           \ 
        |            \
        |             \
        |
        z(depth) negative depths
        
        please notice that okada 92 requiers depth to be defined as negative. Meaning anything under the surface is a negative number. 
        The code takes care of that"""
        
    def GenerateGrid(self,x,y,z):
        
         
         self.xs,self.zs=np.meshgrid(x,z)
         self.gridShape=self.xs.shape
         self.xsFlat = self.xs.flatten()
         self.zsFlat = -1*np.abs(self.zs.flatten())
         self.ysFlat = np.ones_like(self.xsFlat)*self.y
         
#%% trinagluar grid

class gridTriangular(grid):
    def __init__(self,points,cells):
        
        x=points[:,0]
        z=points[:,1]
        
        self.x=ConvertDataToArrayFloat64(x)
        self.z=ConvertDataToArrayFloat64(np.abs(z))
        self.cells=cells
        
        self.xsFlat = self.x
        self.zsFlat = self.z
        self.ysFlat = np.zeros_like(z)
        self.ComputeTriangles()
        
    def ComputeTriangles(self):
        self.triMesh=tri.Triangulation(x=self.xsFlat,y=self.zsFlat,triangles=self.cells)
        
    def Copy(self):
        points=np.zeros([len(self.x),2])
        points[:,0]=self.x
        points[:,1]=self.z
        
        return gridTriangular(points,self.cells)
        
        
#%% function to generate tirg grid

def GenerateTirGridForSimpleDipAngleFault(dipAngle,rectGrid,offsetFromGrid=0.5,mesh_size=0.5,maxExtentFromTrench=None):
    if maxExtentFromTrench is None:
        maxExtentFromTrench=np.max(rectGrid.xs)
        
    faultAtEnd=maxExtentFromTrench*np.tan(np.deg2rad(dipAngle))
    
    
    poly=[[0+offsetFromGrid, 0.0],[maxExtentFromTrench+offsetFromGrid, 0.0],[maxExtentFromTrench+offsetFromGrid,faultAtEnd]]
    
    
    print ("Computing mesh, might take a little bit",flush=True)
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(poly,
            mesh_size=mesh_size)
        mesh = geom.generate_mesh()
        
        
    return gridTriangular(mesh.points,mesh.cells[1].data)

def GenerateTirGridForSimpleDipAngleFaultLowerPlate(dipAngle,rectGrid,offsetFromGrid=0.5,mesh_size=0.5):
    maxExtentFromTrench=np.max(rectGrid.xs)
    faultAtEnd=maxExtentFromTrench*np.tan(np.deg2rad(dipAngle))
    poly=[[0-offsetFromGrid, 0.0],[maxExtentFromTrench-offsetFromGrid, faultAtEnd],[0-offsetFromGrid,faultAtEnd]]
    
    
    print ("Computing mesh, might take a little bit",flush=True)
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(poly,
            mesh_size=mesh_size)
        mesh = geom.generate_mesh()
        
        
    return gridTriangular(mesh.points,mesh.cells[1].data)
    
#%%    
def LoadThurstFaultsAtMidPoint(fautlsDF):
    """ get pandas data frame and return faults object """
    
    return thurstFaultsPositionAtMidFault(xFault=fautlsDF.xFault, yFault=fautlsDF.yFault, zFault=fautlsDF.zFault, faultLengthAlongStrike=fautlsDF.faultLengthAlongStrike, 
                                   faultLengthAlongDip=fautlsDF.faultLengthAlongDip, dipAngle=fautlsDF.dipAngle, strikeAngle=fautlsDF.strikeAngle, dipSlip=fautlsDF.dipSlip)


def LoadNormalFaultsAtMidPoint(fautlsDF):
    """ get pandas data frame and return faults object """
    
    return normalFaultsPositionAtMidFault(xFault=fautlsDF.xFault, yFault=fautlsDF.yFault, zFault=fautlsDF.zFault, faultLengthAlongStrike=fautlsDF.faultLengthAlongStrike, 
                                   faultLengthAlongDip=fautlsDF.faultLengthAlongDip, dipAngle=fautlsDF.dipAngle, strikeAngle=fautlsDF.strikeAngle, dipSlip=fautlsDF.dipSlip)
        
#%%
class thurstFaultsPositionAtMidFault:
  
    """ this class gets location of thurst faults and set the sign of the slip in accordance with the dip angle. It also has a bunch
    of functions to plot the faults. These are the important things to note:
        xFault x location of (west) fault mid point in km
        yFault y location of (north) fault mid point in km
        zFault z location (depth) of fault mid point in km
        faultLengthAlongStrike length of fault along strike in km 
        faultLengthAlongDip length of fault along dip in km 
        dipAngle in deg
        strikeAngle in deg
        dipSlip slip in direction of dip (thurst/normal) in m
        strikeSlip slip in strike slip direction in m
        tensile slip (open/close) in m .
        Lastly it makes sure eveything is in array of float64 because okada92 is very finkey about it"""
        
    def __init__(self,xFault,yFault,zFault, faultLengthAlongStrike,faultLengthAlongDip,dipAngle,strikeAngle,dipSlip,strikeSlip=None,tensileSlip=None): 
        """ this class gets location of thurst faults and set the sign of the slip in accordance with the dip angle. It also has a bunch
        of functions to plot the faults. These are the important things to note:
            xFault x location of (west) fault mid point in km
            yFault y location of (north) fault mid point in km
            zFault z location (depth) of fault mid point in km
            faultLengthAlongStrike length of fault along strike in km 
            faultLengthAlongDip length of fault along dip in km 
            dipAngle in deg
            strikeAngle in deg
            dipSlip slip in direction of dip (thurst/normal) in m
            strikeSlip slip in strike slip direction in m
            tensile slip (open/close) in m .
            Lastly it makes sure eveything is in array of float64 because okada92 is very finkey about it"""
            
   
        self.xFault=ConvertDataToArrayFloat64(xFault)
        self.yFault=ConvertDataToArrayFloat64(yFault)
        self.zFault=ConvertDataToArrayFloat64(zFault)
        self.faultLengthAlongStrike=ConvertDataToArrayFloat64(faultLengthAlongStrike)
        self.faultLengthAlongDip=ConvertDataToArrayFloat64(faultLengthAlongDip)
        self.dipAngle=ConvertDataToArrayFloat64(dipAngle)
        self.strikeAngle=ConvertDataToArrayFloat64(strikeAngle)
            
        dipSlip=self.ChangeSlipDipToMatchFaultType(self.dipAngle,ConvertDataToArrayFloat64(dipSlip))
        self.dipSlip=dipSlip
        self.strikeSlip=self.CheckIfNone(strikeSlip)
        self.tensileSlip=self.CheckIfNone(tensileSlip)
        self.CopyListOfArrays()        
    
            
        self.CheckArrayAreAtTheSameLength(self.listOfallArrays)
            
        self.ShiftPositionOfFaultsToMidPoint()
    
    # def LoadListOfArrays(self,listOfArrays):
        
    #     self.xFault=listOfArrays[0]
    #     self.yFault=listOfArrays[1]
    #     self.zFault=listOfArrays[2]
    #     self.faultLengthAlongStrike=listOfArrays[3]
    #     self.faultLengthAlongDip=listOfArrays[4]
    #     self.dipAngle=listOfArrays[5]
    #     self.strikeAngle=listOfArrays[6]
    #     self.dipSlip=listOfArrays[7]
    #     self.strikeSlip=listOfArrays[8]
    #     self.tensileSlip=listOfArrays[9]
        

    def ReturnDataFrame(self):
        """ reutn pandas dataframe for faults"""
        return pd.DataFrame({'xFault':self.xFault,'yFault':self.yFault,'zFault':self.zFault,'faultLengthAlongStrike':self.faultLengthAlongStrike,
                                      'faultLengthAlongDip':self.faultLengthAlongDip,'dipAngle':self.dipAngle,'strikeAngle':self.strikeAngle,'dipSlip':self.dipSlip,
                                      'strikeSlip':self.strikeSlip,'tensileSlip':self.tensileSlip})
    
    def ComputeMw(self):
        return miscFunctions.GetSlipRutprenLenthDepthReturnMw(self.faultLengthAlongDip*1e3, self.faultLengthAlongStrike*1e3, np.abs(self.dipSlip))
        
    def CopyListOfArrays(self):
        self.listOfallArrays=[self.xFault,self.yFault,self.zFault,self.faultLengthAlongStrike,self.faultLengthAlongDip,
                                           self.dipAngle,self.strikeAngle,
                                           self.dipSlip,self.strikeSlip,self.tensileSlip]
        
    def CheckArrayAreAtTheSameLength(self,listOfArrays):
        for i in range(len(listOfArrays)-1):
            if len(listOfArrays[i+1])!=len(listOfArrays[i]):
                
                raise TypeError("Arrays are not of the same length")
                
            if np.isnan(listOfArrays[i]).any():
                raise TypeError("Arrays contain nan! Not good")
                
            
        if np.isnan(listOfArrays[-1]).any():
            raise TypeError("Arrays contain nan! Not good")
                
                        
        

    def CheckIfNone(self,vector):
        if vector is None:
            vector=np.zeros_like(self.xFault)
        else:
            vector=ConvertDataToArrayFloat64(vector)
            
        return vector
    
    def MergeFaultsStrctureIntoOne(self,listOfFaultsStrctuers):
        
        xFault=np.array([])
        yFault=np.array([])
        zFault=np.array([])
        faultLengthAlongStrike=np.array([])
        faultLengthAlongDip=np.array([])
        dipAngle=np.array([])
        strikeAngle=np.array([])
        dipSlip=np.array([])
        strikeSlip=np.array([])
        tensileSlip=np.array([])
        
        for faultStrcture_i in listOfFaultsStrctuers:
            
            if isinstance(faultStrcture_i,thurstFaultsPositionAtMidFault) is True:
                
                xFault=np.append(xFault,faultStrcture_i.xFault)
                yFault=np.append(yFault,faultStrcture_i.yFault)
                zFault=np.append(zFault,faultStrcture_i.zFault)
                faultLengthAlongStrike=np.append(faultLengthAlongStrike,faultStrcture_i.faultLengthAlongStrike)
                faultLengthAlongDip=np.append(faultLengthAlongDip,faultStrcture_i.faultLengthAlongDip)
                dipAngle=np.append(dipAngle,faultStrcture_i.dipAngle)
                strikeAngle=np.append(strikeAngle,faultStrcture_i.strikeAngle)
                dipSlip=np.append(dipSlip,faultStrcture_i.dipSlip)
                strikeSlip=np.append(strikeSlip,faultStrcture_i.strikeSlip)
                tensileSlip=np.append(tensileSlip,faultStrcture_i.tensileSlip)
                
            else:
                
                raise TypeError("Must be fault type . can't merge otherwise!")
                
        newFault=thurstFaultsPositionAtMidFault(xFault,yFault,zFault,faultLengthAlongStrike,
                                                faultLengthAlongDip,dipAngle,strikeAngle,dipSlip,
                                                strikeSlip,tensileSlip)
            
        return newFault
    
    def SplitFaultAtDipLocation(self,indexOfFault,normalizedPositionToSplit):
        if normalizedPositionToSplit < 0 or normalizedPositionToSplit > 1:
            raise TypeError("This has to be between 0-1 as it is normalizled location")
        
        yFault=self.yFault[indexOfFault]
        faultLengthAlongStrike=self.faultLengthAlongStrike[indexOfFault]
        dipAngle=self.dipAngle[indexOfFault]
        strikeAngle=self.strikeAngle[indexOfFault]
        dipSlip=self.dipSlip[indexOfFault]
        strikeSlip=self.strikeSlip[indexOfFault]
        tensileSlip=self.tensileSlip[indexOfFault]
        
        faultLength=self.faultLengthAlongDip[indexOfFault]
        faultLengthAboveSplittingPoint=faultLength*normalizedPositionToSplit
        faultLengthBelowSplittingPoint=faultLength*normalizedPositionToSplit
        
        x_i,z_i=self.ReturnDeepTipOfault(indexOfFault-1)
        self.DeleteFault(indexOfFault)
        
        faultAboveSplittingPoint=thurstFaultsPositionAtShallowFaultTip(x_i,yFault,z_i,faultLengthAlongStrike,
                                                                       faultLengthAboveSplittingPoint,dipAngle,strikeAngle,dipSlip,strikeSlip,tensileSlip)
                                    
                                              
        x_i,z_i=faultAboveSplittingPoint.ReturnDeepTipOfault()
        faultBelowSplittingPoint=thurstFaultsPositionAtShallowFaultTip(x_i,yFault,z_i,faultLengthAlongStrike,
                                                                       faultLengthBelowSplittingPoint,dipAngle,strikeAngle,dipSlip,strikeSlip,tensileSlip)
        
        newFaults=self.MergeFaultsStrctureIntoOne([faultAboveSplittingPoint,faultBelowSplittingPoint])
        
        return newFaults
                                                  
                                                
        
        
        
    def DeleteFault(self,faultIndex):
        
        self.xFault=np.delete(self.xFault,faultIndex)
        self.yFault=np.delete(self.yFault,faultIndex)
        self.zFault=np.delete(self.zFault,faultIndex)
        self.faultLengthAlongStrike=np.delete(self.faultLengthAlongStrike,faultIndex)
        self.faultLengthAlongDip=np.delete(self.faultLengthAlongDip,faultIndex)
        self.dipAngle=np.delete(self.dipAngle,faultIndex)
        self.strikeAngle=np.delete(self.strikeAngle,faultIndex)
        self.dipSlip=np.delete(self.dipSlip,faultIndex)
        self.strikeSlip=np.delete(self.strikeSlip,faultIndex)
        self.tensileSlip=np.delete(self.tensileSlip,faultIndex)
        
    def ReturnShallowAndDeepTipOfault(self):
        
        for x_i,z_i,L,strikeAngle_i in zip(self.xFault,self.zFault,self.faultLengthAlongDip,self.dipAngle):
            
            dx=L/2*np.cos(np.deg2rad(strikeAngle_i))
            dz=L/2*np.sin(np.deg2rad(strikeAngle_i))
            
            x=[x_i-dx,x_i+dx]
            z=[z_i-dz,z_i+dz]
            
        
        return x,z
    
    def ReturnDeepTipOfault(self,faultIndex=0):
        
        x_i=self.xFault[faultIndex]
        z_i=self.zFault[faultIndex]
        L=self.faultLengthAlongDip[faultIndex]
        strikeAngle_i=self.dipAngle[faultIndex]
     
            
        dx=L/2*np.cos(np.deg2rad(strikeAngle_i))
        dz=L/2*np.sin(np.deg2rad(strikeAngle_i))

             
        return x_i+dx,z_i+dz
            
    def ChangeSlipDipToMatchFaultType(self,dipAngle,slip):

        for i in range(slip.size):
            
            if dipAngle[i] > 0 and dipAngle[i] < 90:
                slip[i]=np.abs(slip[i])
            elif dipAngle[i] >90 and dipAngle[i] <180:
                slip[i]=-1*np.abs(slip[i])
            else:
                raise ValueError ("Dip angle is not defined for the dip " +str(dipAngle))
                
        return slip
        
    def ShiftPositionOfFaultsToMidPoint(self):
        return 
    
    
    def PlotFaultsMapView(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        for x_i,y_i,L,strikeAngle_i in zip(self.xFault,self.yFault,self.faultLengthAlongStrike,self.strikeAngle):
            
            dy=L/2*np.cos(np.deg2rad(strikeAngle_i))
            dx=L/2*np.sin(np.deg2rad(strikeAngle_i))
            
            x=[x_i-dx,x_i+dx]
            y=[y_i-dy,y_i+dy]
            
            ax.plot(x,y)
            
    def PlotFaultsDepthView(self,ax=None,cmap=None,**args):
        if ax is None:
            fig,ax=plt.subplots()
            
        norm = colors.Normalize(vmin=0, vmax=1,)
        
        
        for i,(x_i,z_i,L,strikeAngle_i) in enumerate(zip(self.xFault,self.zFault,self.faultLengthAlongDip,self.dipAngle)):
            
            dx=L/2*np.cos(np.deg2rad(strikeAngle_i))
            dz=L/2*np.sin(np.deg2rad(strikeAngle_i))
            
            x=[x_i-dx,x_i+dx]
            z=[z_i-dz,z_i+dz]
            if cmap is None:
                rgba_color = cm.inferno_r(norm(np.abs(self.dipSlip[i])/np.max(np.abs(self.dipSlip))))
            else:
                rgba_color = cmap(norm(np.abs(self.dipSlip[i])/np.max(np.abs(self.dipSlip))))
                
            cs=ax.plot(x,z,c=rgba_color,**args)
            
        ylim0=ax.get_ylim()[0]
        ylim1=ax.get_ylim()[1]
        
        if ylim1>ylim0:
            ax.invert_yaxis()
            
        return cs
            
    def PlotFaultAtMidPoint(self,ax=None,**args):
        if ax is None:
            fig,ax=plt.subplots()
        
        ax.plot(self.xFault,self.zFault,**args)
        
    def PlotFaultAtDepthAsPoint(self,ax=None,**args):
        if ax is None:
            fig,ax=plt.subplots()
            
        dx=self.faultLengthAlongDip/2*np.cos(np.deg2rad(self.dipAngle))
        dz=self.faultLengthAlongDip/2*np.sin(np.deg2rad(self.dipAngle))
        
        ax.scatter(self.xFault+dz,self.zFault+dz,**args)
            
    def PlotCouplingVsDistance(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            

            
        ax.scatter(self.xFault,1-self.dipSlip/np.max(self.dipSlip),marker='.')
            
        
#%% 
class normalFaultsPositionAtMidFault(thurstFaultsPositionAtMidFault):
    def ChangeSlipDipToMatchFaultType(self,dipAngle,slip):

        for i in range(slip.size):
            
            if dipAngle[i] > 0 and dipAngle[i] < 90:
                slip[i]=-1*np.abs(slip[i])
            elif dipAngle[i] >90 and dipAngle[i] <180:
                slip[i]=np.abs(slip[i])
            else:
                raise ValueError ("Dip angle is not defined for the dip " +str(dipAngle))
                
        return slip
            
            
        
#%%
class thurstFaultsPositionAtShallowFaultTip(thurstFaultsPositionAtMidFault):
    """ this class gets the location of the for the upper tip of the fault and then tranforms
    everything to mid fault so it can work with okada92"""
    
    def ShiftPositionOfFaultsToMidPoint(self):
            
        self.xFault=self.xFault+self.faultLengthAlongDip*np.cos(np.deg2rad(self.dipAngle))/2
        self.zFault=self.zFault+self.faultLengthAlongDip*np.sin(np.deg2rad(self.dipAngle))/2
            
        if self.zFault.any() <0:
            raise ValueError ("Fault is above the surface, change buried depth")
            

class normalFaultsPositionAtShallowFaultTip(thurstFaultsPositionAtShallowFaultTip,normalFaultsPositionAtMidFault):
    pass
            
            
#%%
class thurstFaultsPositionAtDeepFaultTip(thurstFaultsPositionAtMidFault):
    """ this class gets the location of the for the upper tip of the fault and then tranforms
    everything to mid fault so it can work with okada92"""
    
    def ShiftPositionOfFaultsToMidPoint(self):
            
        self.xFault=self.xFault-self.faultLengthAlongDip*np.cos(np.deg2rad(self.dipAngle))/2
        self.zFault=self.zFault-self.faultLengthAlongDip*np.sin(np.deg2rad(self.dipAngle))/2
            
        if self.zFault.any() <0:
            raise ValueError ("Fault is above the surface, change buried depth")
 
class normalFaultsPositionAtDeepFaultTip(thurstFaultsPositionAtDeepFaultTip,normalFaultsPositionAtMidFault):
    pass

#%%
class okadaFor3dFaults:
    """ this class gets grid object and the poistion of the fault in the form of midPointFault object and compute displacmenet.
    for complete doc on how to compute the okada soultion see the doc for function ComputeOkada. Eventuall it returns an xarry for the displacemnt 
    so it is very easy to plot everything. 
    It can also calc stresses with posissionRatio and elastic Modulus that can be set """
    def __init__(self,grid,faultsPositionAtMidPoint,elasticModulus=30.0e9,poissonRatio=0.25,z=None):
        """ grid strcture
        faultsPostion
        elastic moduls
        possion ration
        z"""

        self.faultsPosition=faultsPositionAtMidPoint       
        self.grid=grid
        self.elasticModulus=np.float64(elasticModulus)
        self.poissonRatio=np.float64(poissonRatio)
        self.TestInput(grid)
        self.ComputeOkada()
    
    
    def ComputeOkada(self):
        """    
        # this what the okada92 looks like:
        # u, d, s, flag, flag2 = ok92.okada92(xs, ys, zs, 
        #xc, yc, 
        #depth, length, width, 
        #dip, strike,
        #ss, ds, ts,
        #mu, nu)
        
        # sx ys zs flat vector of all stations location
        #xc and yc are location of the middle of the fault
        #depth = buired depth of the fault at middle of the fault
        #length = length of the whole fault on map view (or trace view)
        # width = fault length aong the the axis of the fault itself
        #dip is angles
        # orenation of the fault
        #ss slip in strike slip
        # ds slip in dip slie (positive is thurst)
        #ts tensile slip
        #mu elastic moduls (usually 30e9) 
        #nu possion ratio
        
        
        #location of the fault is detrmined at middle of the fault. 
        make sure everything is in float64
        """

        xsFlat=self.grid.xsFlat
        ysFlat=self.grid.ysFlat
        zsFlat=self.grid.zsFlat
            
        print("Computing Stresses and displacment",flush=True)
        displacementVector, s, stressVector, flag, flag2 = okada4py.okada92(xsFlat, ysFlat,zsFlat, 
                                            self.faultsPosition.xFault, self.faultsPosition.yFault, self.faultsPosition.zFault, 
                                            self.faultsPosition.faultLengthAlongStrike,  self.faultsPosition.faultLengthAlongDip, 
                                            self.faultsPosition.dipAngle, self.faultsPosition.strikeAngle, 
                                            self.faultsPosition.strikeSlip, self.faultsPosition.dipSlip, self.faultsPosition.tensileSlip,
                                            self.elasticModulus, self.poissonRatio)
        
        #print(self.faultsPosition.dipAngle)
        #
        self.stress=self.ConvertToStresMatrix(stressVector) 
        
        if isinstance(self.grid,elasticDislocations3d.gridxz): ## this is a stupid workaround to not compute displacemnt if I'm using depth view
            
            self.GenerateXarrayFor2D_Displacement(displacementVector,"z",self.grid.z)
        elif isinstance(self.grid,elasticDislocations3d.griddx) or isinstance(self.grid,elasticDislocations3d.grid):
            self.GenerateXarrayFor2D_Displacement(displacementVector,"y",self.grid.y)
            
        elif isinstance(self.grid,elasticDislocations3d.gridTriangular):
            pass
        else:
            print(self.grid)
            raise TypeError("Not sure what to do with this sort of thing!")
            
    def RunForSpeficPoints(self,xsFlat,ysFlat,zsFlat):
        
        xsFlat=ConvertDataToArrayFloat64(xsFlat)
        ysFlat=ConvertDataToArrayFloat64(ysFlat)
        zsFlat=ConvertDataToArrayFloat64(zsFlat)
        
        displacementVector, s, stressVector, flag, flag2 = okada4py.okada92(xsFlat, ysFlat,zsFlat, 
                                            self.faultsPosition.xFault, self.faultsPosition.yFault, self.faultsPosition.zFault, 
                                            self.faultsPosition.faultLengthAlongStrike,  self.faultsPosition.faultLengthAlongDip, 
                                            self.faultsPosition.dipAngle, self.faultsPosition.strikeAngle, 
                                            self.faultsPosition.strikeSlip, self.faultsPosition.dipSlip, self.faultsPosition.tensileSlip,
                                            self.elasticModulus, self.poissonRatio)
        
        return displacementVector, s, stressVector, flag, flag2
        
        
    def ConvertToStresMatrix(self,stressVector):
        stressVector = stressVector.reshape((self.grid.xsFlat.shape[0], 6))
        stress = np.zeros((3, 3, len(self.grid.xsFlat)))
        stress[0,0,:] = stressVector[:,0]
        stress[1,1,:] = stressVector[:,3]
        stress[2,2,:] = stressVector[:,5]
        stress[0,1,:] = stressVector[:,1]
        stress[1,0,:] = stressVector[:,1]
        stress[0,2,:] = stressVector[:,2]
        stress[2,0,:] = stressVector[:,2]
        stress[1,2,:] = stressVector[:,4]
        stress[2,1,:] = stressVector[:,4]
        
        return stress
        
    def GenerateXarrayFor2D_Displacement(self,displacement,extraDimLabel,extraDim):
        
        displacement = displacement.reshape((self.grid.xsFlat.shape[0], 3))
        
        dx=displacement[:,0]
        dy=displacement[:,1]
        dz=displacement[:,2]
        
        dx=dx.reshape(self.grid.gridShape)
        dy=dy.reshape(self.grid.gridShape)
        dz=dz.reshape(self.grid.gridShape)
        
        arrays=[]
        for dd in [dx,dy,dz]:
            array=xr.DataArray(data=dd,dims=["y","x"],coords=[extraDim,self.grid.x])
            array.x.attrs['units']='km'
            array.y.attrs['units']='km'
            array.attrs['units']='m'
            arrays.append(array)
            
        
        self.displacement=xr.Dataset(data_vars={'dz':arrays[2],'dy':arrays[1],'dx':arrays[0]})
        
    def TestInput(self,grid):
        if np.max(grid.zsFlat) > 0 :
            raise ValueError ("Okada92 require under surface values to be negative")
            
    def ReturnStreesesAsObject(self):
        
        if isinstance(grid,elasticDislocations3d.gridxz) is True:
            raise TypeError("This is designed for depth profiles")
            
        Sxx=self.stress[0,0,:].reshape(self.grid.xs.shape)
        Syy=self.stress[2,2,:].reshape(self.grid.xs.shape)
        Sxy=self.stress[0,2,:].reshape(self.grid.xs.shape)
        return ComputeWedge.stresses(Sxx=Sxx, Syy=Syy, X=self.grid.xs*1000, Y=self.grid.zs*1000,
                                    Sxy=Sxy)
    
    def ReturnDisplacemnt(self):
        return self.displacement
        
        
            

        
