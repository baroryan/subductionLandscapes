#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 13:42:02 2023

@author: bar
"""
import multipleFaultsDislocations
import elasticDislocations3d
import numpy as np
import time
import convertContoursToShapes
import gutenbergRichter
import miscClasses
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import miscFunctions
import pandas as pd
from shapely.geometry import Polygon
from shapely.geometry.multipolygon import MultiPolygon
import pygmsh
#%%
def FindColumnbStressChangeAtAngle(Sxx,Syy,Sxy,angle,mu=0.6):

    angle=np.deg2rad(angle)
    c=np.cos(angle)
    s=np.sin(angle)

    return Sxx*s-c*Sxy*(c-mu*s)+(Sxy*s-c*Syy)*(s+mu*c)


def GenerateMesh(dipAngle,offsetFromGrid=1,mesh_size=1):
    gridDepth=elasticDislocations3d.gridxz(x=np.linspace(-200,500.807,5),y=0,z=np.linspace(0,89.460,5))
    mesh=elasticDislocations3d.GenerateTirGridForSimpleDipAngleFault(dipAngle, gridDepth,offsetFromGrid=offsetFromGrid,mesh_size=mesh_size)

    return mesh

def GenerateMeshLowerPlate(dipAngle,offsetFromGrid=1,mesh_size=1):
    gridDepth=elasticDislocations3d.gridxz(x=np.linspace(-200,500.807,5),y=0,z=np.linspace(0,89.460,5))
    mesh=elasticDislocations3d.GenerateTirGridForSimpleDipAngleFaultLowerPlate(dipAngle, gridDepth,offsetFromGrid=offsetFromGrid,mesh_size=mesh_size)

    return mesh


def ComputeNormalizedColumbStressChange(dipAngle,lockingDistanceAlongFault,LockingWidthAlongFault,mesh,alpha=0.2,mu=0.6,interfaceFaults=None):
    if interfaceFaults is None:
        linearSlipZeroToOneOver=multipleFaultsDislocations.linearSlipDistribution(numberOfFaultToSlip=100,maxNormalizedSlip=1,minNormalizedSlip=0)
        chileMegathurstFault=multipleFaultsDislocations.subsequentDislocationsToBuildComplexGeometries(dipAngles=[dipAngle,dipAngle,dipAngle],lengthsOfDislocations=[lockingDistanceAlongFault,LockingWidthAlongFault,50000],depthOfFirstDislocationShallowTip=0,distanceOfFirstDislocationShallowTip=0,
                                                                                                       slipDistributionObjects=[None,linearSlipZeroToOneOver,None],dipSlips=[0,1e-3,1e-3])
    
        chileMegathurstFault=chileMegathurstFault.splittedDislocations
        
    else:
        chileMegathurstFault=interfaceFaults
        
    gridDepth=elasticDislocations3d.gridxz(x=np.linspace(-200,500.807,5),y=0,z=np.linspace(0,89.460,5))
    chileMegathurstFaultObject=elasticDislocations3d.okadaFor3dFaults(gridDepth,chileMegathurstFault)


    _, _, stressVector, _, _=chileMegathurstFaultObject.RunForSpeficPoints(np.array(mesh.xsFlat,dtype=np.float64),np.zeros_like(mesh.xsFlat),np.array(-1*mesh.zsFlat,dtype=np.float64))

    stress=stressVector.reshape([len(mesh.xsFlat),6])
    coulmbStressChange=FindColumnbStressChangeAtAngle(Sxx=stress[:,0],Syy=stress[:,5],Sxy=stress[:,2],angle=30,mu=mu-alpha)
    coulmbStressChange=np.where(coulmbStressChange>0,0,coulmbStressChange)
    coulmbStressChange=np.abs(coulmbStressChange)
    coulmbStressChange=(coulmbStressChange/np.nanmax(coulmbStressChange))
    
    
    return coulmbStressChange


def ComputeColumbStressChange(dipAngle,lockingDistanceAlongFault,LockingWidthAlongFault,mesh,alpha=0.2,mu=0.6,interfaceFaults=None,compressive=True):
    if interfaceFaults is None:
        linearSlipZeroToOneOver=multipleFaultsDislocations.linearSlipDistribution(numberOfFaultToSlip=100,maxNormalizedSlip=1,minNormalizedSlip=0)
        chileMegathurstFault=multipleFaultsDislocations.subsequentDislocationsToBuildComplexGeometries(dipAngles=[dipAngle,dipAngle,dipAngle],lengthsOfDislocations=[lockingDistanceAlongFault,LockingWidthAlongFault,50000],depthOfFirstDislocationShallowTip=0,distanceOfFirstDislocationShallowTip=0,
                                                                                                       slipDistributionObjects=[None,linearSlipZeroToOneOver,None],dipSlips=[0,1e-3,1e-3])
    
        chileMegathurstFault=chileMegathurstFault.splittedDislocations
        
    else:
        chileMegathurstFault=interfaceFaults
        
    gridDepth=elasticDislocations3d.gridxz(x=np.linspace(-200,500.807,5),y=0,z=np.linspace(0,89.460,5))
    chileMegathurstFaultObject=elasticDislocations3d.okadaFor3dFaults(gridDepth,chileMegathurstFault)


    _, _, stressVector, _, _=chileMegathurstFaultObject.RunForSpeficPoints(np.array(mesh.xsFlat,dtype=np.float64),np.zeros_like(mesh.xsFlat),np.array(-1*mesh.zsFlat,dtype=np.float64))

    stress=stressVector.reshape([len(mesh.xsFlat),6])
    
    
    if compressive is True:
        coulmbStressChange=FindColumnbStressChangeAtAngle(Sxx=stress[:,0],Syy=stress[:,5],Sxy=stress[:,2],angle=30,mu=mu-alpha)
        coulmbStressChange=np.where(coulmbStressChange>0,0,coulmbStressChange)
    elif compressive is False:
        coulmbStressChange=FindColumnbStressChangeAtAngle(Sxx=stress[:,0],Syy=stress[:,5],Sxy=stress[:,2],angle=60,mu=mu-alpha)
        coulmbStressChange=np.where(coulmbStressChange<0,0,coulmbStressChange)
    
    
    return coulmbStressChange


def ComputeSurfaceDisplacement(dipAngle,lockingDistanceAlongFault,LockingWidthAlongFault):
    
    linearSlipZeroToOneOver=multipleFaultsDislocations.linearSlipDistribution(numberOfFaultToSlip=100,maxNormalizedSlip=1,minNormalizedSlip=0)
    chileMegathurstFault=multipleFaultsDislocations.subsequentDislocationsToBuildComplexGeometries(dipAngles=[dipAngle,dipAngle,dipAngle],lengthsOfDislocations=[lockingDistanceAlongFault,LockingWidthAlongFault,50000],depthOfFirstDislocationShallowTip=0,distanceOfFirstDislocationShallowTip=0,
                                                                                                   slipDistributionObjects=[None,linearSlipZeroToOneOver,None],dipSlips=[0,1e-3,1e-3])

    chileMegathurstFault=chileMegathurstFault.splittedDislocations
    grid=elasticDislocations3d.grid(x=np.linspace(-200,500.807,1000), y=np.array([0]))
    chileMegathurstFaultObject=elasticDislocations3d.okadaFor3dFaults(grid,chileMegathurstFault)


    return chileMegathurstFaultObject.ReturnDisplacemnt()

def ReturnInterfaceFaults(dipAngle,lockingDistanceAlongFault,LockingWidthAlongFault):
    
    linearSlipZeroToOneOver=multipleFaultsDislocations.linearSlipDistribution(numberOfFaultToSlip=100,maxNormalizedSlip=1,minNormalizedSlip=0)
    chileMegathurstFault=multipleFaultsDislocations.subsequentDislocationsToBuildComplexGeometries(dipAngles=[dipAngle,dipAngle,dipAngle],lengthsOfDislocations=[lockingDistanceAlongFault,LockingWidthAlongFault,50000],depthOfFirstDislocationShallowTip=0,distanceOfFirstDislocationShallowTip=0,
                                                                                                   slipDistributionObjects=[None,linearSlipZeroToOneOver,None],dipSlips=[0,1e-3,1e-3])

    chileMegathurstFault=chileMegathurstFault.splittedDislocations
    
    return chileMegathurstFault


def GenerateRandomFloatBetweenMinAndMax(minLim,maxLim,N=1):
        return minLim+(maxLim-minLim)*np.random.rand(N)

def ConvertDataToContuors(data,contours):
    contours=np.asarray(contours)
    for i in range(len(contours)-1):
        data=np.where((data >contours[i] ) & (data <contours[i+1] ),contours[i],data)

    return data

def ReturnSlip(Mw):
    """ Get Mw returns slip, fig 11 Wells & coppermith et al., 1994"""
    return 10**((Mw-6.93)/0.82)


def RetrunMwForDepthExtent(D):
    """Get depth exten in km returns Mw. figure 15 Wells & Coppersmith et al., 1994"""
    return 4.06+2.25*np.log10(D)

def ReturnDepthExtent(Mw):
    """Get Mw return depth exten in km . figure 15 Wells & Coppersmith et al., 1994"""
    return 10**((Mw-4.06)/2.25)

def ReturnRuptureLength(Mw):
    """Get Mw return surface rupture length  in km . figure 9 Wells & Coppersmith et al., 1994"""
    return 10**((Mw-4.38)/1.49)

def GenerateBinodalRandInt(oneInt,secondInt,N):
    rand=np.random.rand(N)
    rand=np.where(rand>0.5,int(oneInt),int(secondInt))

    return rand

def ComputeSlip(Mw,area,elasticModuls=30e9):
    """ get Mw area in sq m and return slip in meters"""
    return 10**((Mw*1.5+9.049))/(30e9*area)

def ComputeAlongStrikeRupture(Mw,depthExtent,slip,elasticModuls=30e9):
    """ get Mw slip and depth extent and return along strike extent"""
    return 10**((Mw*1.5+9.049))/(elasticModuls*slip*depthExtent)

#...............................................................................
def sumtriangles( xy, z, triangles ):
    """ integrate scattered data, given a triangulation
    zsum, areasum = sumtriangles( xy, z, triangles )
    In:
        xy: npt, dim data points in 2d, 3d ...
        z: npt data values at the points, scalars or vectors
        triangles: ntri, dim+1 indices of triangles or simplexes, as from
http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html
    Out:
        zsum: sum over all triangles of (area * z at midpoint).
            Thus z at a point where 5 triangles meet
            enters the sum 5 times, each weighted by that triangle's area / 3.
        areasum: the area or volume of the convex hull of the data points.
            For points over the unit square, zsum outside the hull is 0,
            so zsum / areasum would compensate for that.
            Or, make sure that the corners of the square or cube are in xy.
    """
        # z concave or convex => under or overestimates
    npt, dim = xy.shape
    ntri, dim1 = triangles.shape
    assert npt == len(z), "shape mismatch: xy %s z %s" % (xy.shape, z.shape)
    assert dim1 == dim+1, "triangles ? %s" % triangles.shape
    zsum = np.zeros( z[0].shape )
    areasum = 0
    dimfac = np.prod( np.arange( 1, dim+1 ))
    for tri in triangles:
        corners = xy[tri]
        t = corners[1:] - corners[0]
        if dim == 2:
            area = abs( t[0,0] * t[1,1] - t[0,1] * t[1,0] ) / 2
        else:
            area = abs( np.linalg.det( t )) / dimfac  # v slow
        zsum += area * z[tri].mean(axis=0)
        areasum += area
    return (zsum, areasum)

def GetMeshColumbReturnSumOfTriangels(mesh,columnb):
    meshPoints=np.zeros([len(mesh.xsFlat),2])
    meshPoints[:,0]=mesh.triMesh.x
    meshPoints[:,1]=mesh.triMesh.y
    sumCol,_=sumtriangles(meshPoints,columnb,mesh.triMesh.triangles)

    return sumCol

def ComputeInterpolation(mesh,columnbStressChange,kindOfInterpolation='geom'):
    
    meshPoints=np.zeros([len(mesh.xsFlat),2])
    meshPoints[:,0]=mesh.triMesh.x
    meshPoints[:,1]=mesh.triMesh.y
    sumCol=sumtriangles(meshPoints,columnbStressChange,mesh.triMesh.triangles)
    interp_cubic_geom = mtri.CubicTriInterpolator(mesh.triMesh, columnbStressChange/sumCol[0], kind=kindOfInterpolation)
    
    return interp_cubic_geom

def GetCutoffReturnContourPolygon(mesh,upperPlateColumnb,cutoff,index=None):
    """ this function computes the cutoff polygon based on cutoff (which is a scalar value)
    for some reason for the simple cases the convert polygon to shaeply package with a simple wedge index = -1 works
    but the Himalayas case where I have a complicated wedge I need to use index=0.
    Update now this functions is looking for the polygon with the largest area and uses it as the one to use. unless I give it specficall
    updae and hopefully this is the last one. May 3rd 2023
    - now I sort values based on columnb stress change and run a while loop so the first point that returns a polygon is that one I use """
    
    cutoffContours=np.where(upperPlateColumnb>cutoff,100,0)
    cs=plt.tricontourf(mesh.triMesh,cutoffContours,levels=2)
    plt.close()
    cutOffContour=convertContoursToShapes.getContourGenerateSheplyContours(cs)
    
    if index is None:
        
        indexesOfPoints=np.argsort(np.abs(upperPlateColumnb-0.5))
        columnbIndex=0
        poly=cutOffContour.FindPointInWhatContourAndReturnPolygon(mesh.xsFlat[indexesOfPoints[columnbIndex]],mesh.zsFlat[indexesOfPoints[columnbIndex]])

        while isinstance(poly, Polygon) is False and isinstance(poly, MultiPolygon) is False:
            columnbIndex=columnbIndex+1
            poly=cutOffContour.FindPointInWhatContourAndReturnPolygon(mesh.xsFlat[indexesOfPoints[columnbIndex]],mesh.zsFlat[indexesOfPoints[columnbIndex]])
            
        cutoffPolygon=cutOffContour.ConvertPolyExteriorPoints(poly)
        
    else:
        index=int(index)
        cutoffPolygon=cutOffContour.ConvertPolyExteriorPoints(cutOffContour.contours[index].poly)
    

    cutoffPolygon=cutoffPolygon.get_xy()
    
    
    return cutoffPolygon
    

def SetRupturePolygon(cutoffPolygon):
    
    index=np.argmax(cutoffPolygon[:,1])
    xForDeepestPoint=cutoffPolygon[index,0]
    bigPoly=miscClasses.wedgeGeo(theta=30,h=np.max(cutoffPolygon[:,1]))
    
    rupturePoly=np.array([[0,0],[float(bigPoly.x+xForDeepestPoint),0],[xForDeepestPoint,np.max(cutoffPolygon[:,1])]])
    
    return rupturePoly
    

#%% function releated to population EQs

def CheckIfTooShallowToRptureOld(depth):
    cutoffDepth=2
    chanceToPass=0.05
    rejectionNumber=np.random.rand(len(depth))
    prefectOfCutoff=(1-(cutoffDepth-depth)/cutoffDepth)*(2*chanceToPass)
    prefectOfCutoff=np.where(depth>2,10,prefectOfCutoff)
    rejection=np.where( prefectOfCutoff>rejectionNumber,True,False) 
    
    return rejection

def CheckIfTooShallowToRpture(depth):
    """ this functions returns only 5% of those values between cutoffDepth and zero. and give up on values under zero"""    
    cutoffDepth=2
    allDepth=pd.DataFrame({'depth':depth})
    inRangeDepth=allDepth.loc[ (allDepth.depth < cutoffDepth) & (allDepth.depth > 0) ]
    indexToBeIncluded=inRangeDepth.sample(int(len(inRangeDepth)*0.05)).index
    mask=np.zeros_like(depth,dtype=bool)
    mask=np.where(depth>2,True,mask)
    mask[indexToBeIncluded]=True

    return mask


def ReturnShallowPoints(xdeep,ydeep,angle,L):

    xShallow=xdeep-np.sin(np.deg2rad(90-angle))*L # converting to shallow points of the fault
    yShallow=ydeep-np.abs(np.cos(np.deg2rad(90-angle))*L)
    #shallow=yShallow>2

    #return xShallow,yShallow,shallow
    return xShallow,yShallow

def ReturnDeepPoints(xshallow,yshallow,angle,L):
    xDeep=xshallow-np.sin(np.deg2rad(angle))*L # converting to shallow points of the fault
    yDeep=yshallow+np.abs(np.cos(np.deg2rad(angle))*L)
    #shallow=yrand>2

    #return xDeep,yDeep,shallow
    return xDeep,yDeep

#%%
class ReturnCompressiveEarthquackePositionAngleWithShallow:
    def __init__(self):
        pass
    
    def ReturnRandomAngle(self,N):
        return GenerateBinodalRandInt(30,150,N)
    
    def CheckIfTooShallowToRpture(self,depth):
        return CheckIfTooShallowToRpture(depth)
    
    def GenerateUnifromRandomEQs(self,rupturePolygon,N):
        xmin=np.min(rupturePolygon[:,0]);xmax=np.max(rupturePolygon[:,0]);ymin=np.min(rupturePolygon[:,1]);ymax=np.max(rupturePolygon[:,1])
        xrand=GenerateRandomFloatBetweenMinAndMax(xmin,xmax,N)
        yrand=GenerateRandomFloatBetweenMinAndMax(ymin,ymax,N)
        inRupturePolygon,_=miscFunctions.CheckIfXandYptsInsidePolygon(xrand, yrand, rupturePolygon)
    
        xrand=xrand[inRupturePolygon];yrand=yrand[inRupturePolygon]
        
        return xrand,yrand
        
    
    def ReturnEarthQuckesPositionAndAngle(self,Mw,rupturePolygon,wedgePolygon,interpolationForChance):
        
        xrand,yrand=self.GenerateUnifromRandomEQs(rupturePolygon, len(Mw))
        N=len(xrand)
        Mw=Mw[0:N]

        depthExtentOfEQs=ReturnDepthExtent(Mw)
        #angleRandForComputingShallowPoint=90-(np.abs(GenerateBinodalRandInt(0,180,N)-angleInTheRegion))
        angleRandForComputingShallowPoint=self.ReturnRandomAngle(N)
    
    
        #xShallow,yShallow,shallow=CheckIfShallow(xrand,yrand,angleRandForComputingShallowPoint,depthExtent)
        xShallow,yShallow=ReturnShallowPoints(xrand,yrand,angleRandForComputingShallowPoint,depthExtentOfEQs)
        shallowToRpture=self.CheckIfTooShallowToRpture(yShallow)
    
        angleRandForComputingShallowPoint=angleRandForComputingShallowPoint[shallowToRpture];xrand=xrand[shallowToRpture];yrand=yrand[shallowToRpture];depthExtentOfEQs=depthExtentOfEQs[shallowToRpture];Mw=Mw[shallowToRpture];xShallow=xShallow[shallowToRpture];yShallow=yShallow[shallowToRpture]
    
    
    
    
        shallowisin,_=miscFunctions.CheckIfXandYptsInsidePolygon(xShallow,yShallow, wedgePolygon)
    
        angleRandForComputingShallowPoint=angleRandForComputingShallowPoint[shallowisin];xrand=xrand[shallowisin];yrand=yrand[shallowisin];depthExtentOfEQs=depthExtentOfEQs[shallowisin];Mw=Mw[shallowisin]
    
        try :
            chance=interpolationForChance(xrand,yrand)
            
        except ValueError :
            arrayForRecInterpolation=np.zeros([len(xrand),2])
            arrayForRecInterpolation[:,0]=xrand;arrayForRecInterpolation[:,1]=yrand
            chance=interpolationForChance(arrayForRecInterpolation)
    
        rejectionNumber=np.random.rand(len(xrand))
    
        accepteanceIndex=np.where(chance.data>rejectionNumber,True,False)



        return xrand[accepteanceIndex],yrand[accepteanceIndex],angleRandForComputingShallowPoint[accepteanceIndex],Mw[accepteanceIndex],depthExtentOfEQs[accepteanceIndex]
    
#%%
class ReturnCompressiveEarthquackePositionAngleWithShallow_40(ReturnCompressiveEarthquackePositionAngleWithShallow):
    def ReturnRandomAngle(self,N):
        return GenerateBinodalRandInt(30+10,150-10,N)
    
class ReturnCompressiveEarthquackePositionAngleWithShallow_20(ReturnCompressiveEarthquackePositionAngleWithShallow):
    def ReturnRandomAngle(self,N):
        return GenerateBinodalRandInt(30-10,150+10,N)
    
class ReturnCompressiveEarthquackePositionAngleWithShallow_std(ReturnCompressiveEarthquackePositionAngleWithShallow):
    def ReturnRandomAngle(self,N):
        std = np.random.normal(0, 10, N)
        
        return GenerateBinodalRandInt(30,150,N)+std
    
#%%
class ReturnCompressiveEarthquackePositionAngleWithNoShallow(ReturnCompressiveEarthquackePositionAngleWithShallow):
    def CheckIfTooShallowToRpture(self,depth):
        return np.ones_like(depth,dtype=bool)

#%%  
class ReturnExtensionalEarthquackePositionAngle(ReturnCompressiveEarthquackePositionAngleWithShallow):
    def ReturnRandomAngle(self,N):
        return GenerateBinodalRandInt(60,120,N)
#%%    

def ReturnEarthQuckesPositionAndAngleWithRejection(Mw,rupturePolygon,wedgePolygon,interpolationForChance):
    xmin=np.min(rupturePolygon[:,0]);xmax=np.max(rupturePolygon[:,0]);ymin=np.min(rupturePolygon[:,1]);ymax=np.max(rupturePolygon[:,1])

    MwRejected=[]
    
    N=len(Mw)
    xrand=GenerateRandomFloatBetweenMinAndMax(xmin,xmax,N)
    yrand=GenerateRandomFloatBetweenMinAndMax(ymin,ymax,N)

    inRupturePolygon,_=miscFunctions.CheckIfXandYptsInsidePolygon(xrand, yrand, rupturePolygon)
    MwRejected.append(Mw[~inRupturePolygon])

    xrand=xrand[inRupturePolygon];yrand=yrand[inRupturePolygon];N=len(xrand);Mw=Mw[inRupturePolygon]




    depthExtentOfEQs=ReturnDepthExtent(Mw)
    #angleRandForComputingShallowPoint=90-(np.abs(GenerateBinodalRandInt(0,180,N)-angleInTheRegion))
    angleRandForComputingShallowPoint=GenerateBinodalRandInt(30,150,N)


    #xShallow,yShallow,shallow=CheckIfShallow(xrand,yrand,angleRandForComputingShallowPoint,depthExtent)
    xShallow,yShallow=ReturnShallowPoints(xrand,yrand,angleRandForComputingShallowPoint,depthExtentOfEQs)
    shallowToRpture=CheckIfTooShallowToRpture(yShallow)

    angleRandForComputingShallowPoint=angleRandForComputingShallowPoint[shallowToRpture];xrand=xrand[shallowToRpture];yrand=yrand[shallowToRpture];depthExtentOfEQs=depthExtentOfEQs[shallowToRpture];Mw=Mw[shallowToRpture];xShallow=xShallow[shallowToRpture];yShallow=yShallow[shallowToRpture]
    MwRejected.append(Mw[~shallowToRpture])



    shallowisin,_=miscFunctions.CheckIfXandYptsInsidePolygon(xShallow,yShallow, wedgePolygon)
    MwRejected.append(Mw[~shallowisin])
    
    angleRandForComputingShallowPoint=angleRandForComputingShallowPoint[shallowisin];xrand=xrand[shallowisin];yrand=yrand[shallowisin];depthExtentOfEQs=depthExtentOfEQs[shallowisin];Mw=Mw[shallowisin]


    chance=interpolationForChance(xrand,yrand)

    rejectionNumber=np.random.rand(len(xrand))

    accepteanceIndex=np.where(chance.data>rejectionNumber,True,False)
    MwRejected.append(Mw[~accepteanceIndex])


    return xrand[accepteanceIndex],yrand[accepteanceIndex],angleRandForComputingShallowPoint[accepteanceIndex],Mw[accepteanceIndex],depthExtentOfEQs[accepteanceIndex],np.array(MwRejected)


def CheckStdToMeanRation(dz):
    """ check the ration of std/mean for max std point"""
    
    std=dz.std('y').values
    mean=dz.mean('y').values
    
    index=np.argmax(std)
    return std[index]/np.abs(mean[index])
    
     


def ComputeSlipAndLength(Mw,depthExtent):
    faultLenthAlongStrike=ReturnRuptureLength(Mw)
    slip=ComputeSlip(Mw,faultLenthAlongStrike*depthExtent*1e6)
    
    return slip,faultLenthAlongStrike



def ComputeMeshAndIntepolation(dipAngle,distanceLocked,distanceSemiLocked,cutoff,mesh=None,offsetFromGrid=1,mesh_size=1,alpha=0.2,mu=0.6):

        
    if mesh is None:
        upperPlateMesh=GenerateMesh(dipAngle,offsetFromGrid=offsetFromGrid,mesh_size=mesh_size)
    else:
        upperPlateMesh=mesh
            
    smallWedge=miscClasses.wedgeGeo(theta=dipAngle,x=distanceLocked)
    bigWedge=miscClasses.wedgeGeo(theta=dipAngle,x=distanceLocked+distanceSemiLocked)
    
    upperPlateColumnb=ComputeNormalizedColumbStressChange(dipAngle,smallWedge.R,bigWedge.R-smallWedge.R,upperPlateMesh,alpha=alpha,mu=mu)
    upperPlateColumnb=np.where(upperPlateColumnb>cutoff,upperPlateColumnb,0)
    interp_cubic_geom=ComputeInterpolation(upperPlateMesh,upperPlateColumnb)
    cutoffPolygon=GetCutoffReturnContourPolygon(upperPlateMesh,upperPlateColumnb,cutoff)
        
    return cutoffPolygonAndInterpolation(cutoffPolygon=cutoffPolygon,intepolation=interp_cubic_geom,upperPlateColumnb=upperPlateColumnb,mesh=upperPlateMesh)
        
#%%    
class cutoffPolygonAndInterpolation:
    def __init__(self,cutoffPolygon,intepolation,upperPlateColumnb=None,mesh=None):
        self.cutoffPolygon=cutoffPolygon
        self.intepolation=intepolation
        self.upperPlateColumnb=upperPlateColumnb
        self.mesh=mesh
        
    def ComputeInterpolation(self):
        self.intepolation=ComputeInterpolation(self.mesh,self.upperPlateColumnb)
        
    def ReturnObjects(self):
        return self.intepolation,self.cutoffPolygon
#%%
class guessEQ:
        
    def GuessEQs(self,cutoffAndIntpolationObject,maxMw,lengthOfDomain):
        
        interp_cubic_geom,cutoffPolygon=cutoffAndIntpolationObject.ReturnObjects()
        Mw=gutenbergRichter.RandomSamplerForGR(b=self.b, minMw=self.minEQsize, maxMw=maxMw,sampleSize=self.sampleSizeForOneRun)
        
        poistionOfEQ=self.CreateTypeOfEQs()
        xDeep,yDeep,angle,Mw,depthExtent=poistionOfEQ.ReturnEarthQuckesPositionAndAngle(Mw,cutoffPolygon,cutoffPolygon,interp_cubic_geom)
        slip,alongStrikeLength=ComputeSlipAndLength(Mw,depthExtent)
        
        
        zPosition=GenerateRandomFloatBetweenMinAndMax(0,lengthOfDomain,len(angle))


        faults=elasticDislocations3d.thurstFaultsPositionAtDeepFaultTip(xDeep,zPosition,yDeep,alongStrikeLength,depthExtent,angle,np.zeros_like(angle),slip)
        faultsDataFrame=faults.ReturnDataFrame()
        
        return faultsDataFrame
    
    def CreateTypeOfEQs(self):
        return ReturnCompressiveEarthquackePositionAngleWithShallow()
    
    def GuessUniformRandomEQs(self,cutoffAndIntpolationObject,N):
        _,cutoffPolygon=cutoffAndIntpolationObject.ReturnObjects()
        poistionOfEQ=self.CreateTypeOfEQs()
        xrand,yrand=poistionOfEQ.GenerateUnifromRandomEQs(cutoffPolygon,N)
        
        return xrand,yrand
        
#%%
class guessEQNoShallow(guessEQ):
    def CreateTypeOfEQs(self):
        return ReturnCompressiveEarthquackePositionAngleWithNoShallow()
    

#%%
class longTermUplift(guessEQ):
    def __init__(self,sampleSizeForOneRun=1e7,minRationOfMeanToStd=0.2,a=1,b=0.9,maxIteration=100,maxSizeOfDomain=5,dx=4,ny=70,minEQsize=5,maxEQsize=None,maxNumOfEQs=None):

        self.sampleSizeForOneRun=sampleSizeForOneRun
        self.minRationOfMeanToStd=minRationOfMeanToStd
        self.a=a
        self.b=b
        self.maxIteration=maxIteration
        #self.sampleSizeForOneRun=1e7
        self.maxSizeOfDomain=maxSizeOfDomain
        self.dx=dx
        self.ny=ny
        self.minEQsize=minEQsize
        self.maxNumOfEQs=maxNumOfEQs
        self.maxEQsize=maxEQsize
        
        
        #finalDisplacment,allFaultsDataFrame,ratioOfMeanToStd=self.ComputeLongTermDisplacement(cutoffAndIntpolationObject)
        
        #self.finalDisplacment=finalDisplacment
        #self.allFaultsDataFrame=allFaultsDataFrame
        #self.ratioOfMeanToStd=ratioOfMeanToStd
        
    # def ReturnDisplacemnt(self):
    #     return self.finalDisplacment
    
    def SetDomainLength(self,maxMw):
        return ReturnRuptureLength(maxMw)*self.maxSizeOfDomain

    
    def GetMaxMw(self,cutoffAndIntpolationObject):
        
        if self.maxEQsize is not None:
            return self.maxEQsize
        
        interp_cubic_geom,cutoffPolygon=cutoffAndIntpolationObject.ReturnObjects()
        reallyBigWedge=miscClasses.wedgeGeo(theta=30,h=np.max(cutoffPolygon[:,1]))
        maxMw=RetrunMwForDepthExtent(reallyBigWedge.R)
        
        return maxMw
    
    def GetGrid(self,lengthOfDomain):
        grid=elasticDislocations3d.grid(x=np.arange(0,253,self.dx),y=np.linspace(0,lengthOfDomain,self.ny))
        return grid
    
    
    def LoadFaults(self,faultsDataFrame):
        faults=elasticDislocations3d.LoadThurstFaultsAtMidPoint(faultsDataFrame)
        return faults
        
    def ComputeLongTermDisplacement(self,cutoffAndIntpolationObject):
        """ compute long term uplift for dipAngle,distanceLocked,distaneSemi locked and cutoff, can also use a mesh given from above.
        return xarrays object. Can also give sample size
        This fuctions runs until cutoff is reached. But it does it compution on very low res array and then compute it on high res array for all EQs"""
    
        t1=time.time()

        maxMw=self.GetMaxMw(cutoffAndIntpolationObject)
        lengthOfDomain=self.SetDomainLength(maxMw)
        grid=self.GetGrid(lengthOfDomain)
        totalDisplacment=grid.GenerateEmptyXRGrid()
        ratioOfMeanToStd=1000
        
        c=0
        faultsList=[]
        totalNumOfEQs=0
        
        while ratioOfMeanToStd>self.minRationOfMeanToStd and c<self.maxIteration:
            print("adding more EQs because ratioOfMeanToStd is: "+str(np.round(ratioOfMeanToStd*100,2))+" and I only have "+str(totalNumOfEQs)+ " EQs",flush='True')
            
            faultsDataFrame=self.GuessEQs(cutoffAndIntpolationObject, maxMw, lengthOfDomain)
            totalNumOfEQs=totalNumOfEQs+len(faultsDataFrame)
            
            if self.maxNumOfEQs is not None and self.maxNumOfEQs>0: 
                """ this part just make sure I provide a set number of EQs and not above"""
                if totalNumOfEQs>self.maxNumOfEQs:
                    index=self.maxNumOfEQs-(totalNumOfEQs-len(faultsDataFrame))-1
                    print ("Reached max num of EQs, this is the last iteration, only using "+ str(index)+ " EQs out of "+ str(len(faultsDataFrame))+ " avilable at this iteration",flush=True)
                    faultsDataFrame=faultsDataFrame.loc[0:index].copy()
                    c=self.maxIteration+10
                
                
            faultsDataFrame['counter']=c
            
            faultsList.append(faultsDataFrame)
            
            #faults=elasticDislocations3d.LoadThurstFaultsAtMidPoint(faultsDataFrame)
            faults=self.LoadFaults(faultsDataFrame)
            deformation=elasticDislocations3d.okadaFor3dFaults(grid,faults)
            displacment=deformation.ReturnDisplacemnt()
            
            totalDisplacment=displacment+totalDisplacment
            
            ratioOfMeanToStd=CheckStdToMeanRation(totalDisplacment.dz)
            
            c=c+1
        
       
        allFaultsDataFrame=pd.concat(faultsList,ignore_index=True).reset_index(drop=True)

        
        print("Finsihed compuation in mins: " + str(np.floor((time.time()-t1)/60)),flush='True')
        
        return totalDisplacment,allFaultsDataFrame,ratioOfMeanToStd
#%%
class guessEQsOnly(longTermUplift,guessEQNoShallow):
    def GuessNumOfEQ(self,cutoffAndIntpolationObject,N=None):
        
        if N is None:
            maxNumOfEQs=self.maxNumOfEQs
        else:
            maxNumOfEQs=N
        
    
        maxMw=self.GetMaxMw(cutoffAndIntpolationObject)
        lengthOfDomain=self.SetDomainLength(maxMw)
        
        c=0
        faultsDataFrame=[] 
        faultsDataFrameList=[]
        

            
        while len(faultsDataFrame)<maxNumOfEQs and c<self.maxIteration:
            faultsDataFrameList.append(self.GuessEQs(cutoffAndIntpolationObject, maxMw, lengthOfDomain))
            faultsDataFrame=self.GetListOfPandasReturnPanda(faultsDataFrameList)
            c=c+1
            
           
            
        if faultsDataFrame.empty:
            raise ValueError ("couldn't guess enough EQ")
            
            
        try :    
            return faultsDataFrame.loc[0:maxNumOfEQs-1]
        except:
            return faultsDataFrame
    
    def GuessRandomUniformEQs(self,cutoffAndIntpolationObject,N=None):
        if N is not None:
            maxNumOfEQs=int(N)
        else:
            maxNumOfEQs=int(self.maxNumOfEQs)
        
        xrand=np.array([])
        yrand=np.array([])
        c=0
        while len(xrand)<maxNumOfEQs and c<self.maxIteration:
            
            x,y=self.GuessUniformRandomEQs(cutoffAndIntpolationObject,2*self.maxNumOfEQs)
            xrand=np.append(xrand, x)
            yrand=np.append(yrand, y)
            c=c+1
            
        if c>self.maxIteration:
            raise ValueError ("couldn't guess enough EQ")
        
        position=np.zeros([len(xrand),2])
        position[:,0]=xrand;position[:,1]=yrand
        print(maxNumOfEQs)
        return position[0:maxNumOfEQs]

    def GetListOfPandasReturnPanda(self,listOfPandas):    
        return pd.concat(listOfPandas,ignore_index=True)
    
    def ReturnDistanceAndDepth(self,faultsDataFrame):
        return faultsDataFrame.loc[:,['xFault','zFault']].values
        
        
        
    
#%%

class longTermUpliftFixedLengthOfDomain(longTermUplift):
    def SetDomainLength(self,maxMw):
        return self.maxSizeOfDomain
    
    
#%%
class longTermUpliftExtensional(longTermUplift):
    def GuessEQs(self,cutoffAndIntpolationObject,maxMw,lengthOfDomain):
        
        interp_cubic_geom,cutoffPolygon=cutoffAndIntpolationObject.ReturnObjects()
        Mw=gutenbergRichter.RandomSamplerForGR(b=self.b, minMw=self.minEQsize, maxMw=maxMw,sampleSize=self.sampleSizeForOneRun)
        poistionOfEQ=ReturnExtensionalEarthquackePositionAngle()
        xDeep,yDeep,angle,Mw,depthExtent=poistionOfEQ.ReturnEarthQuckesPositionAndAngle(Mw,cutoffPolygon,cutoffPolygon,interp_cubic_geom)
        slip,alongStrikeLength=ComputeSlipAndLength(Mw,depthExtent)
        
        
        zPosition=GenerateRandomFloatBetweenMinAndMax(0,lengthOfDomain,len(angle))


        faults=elasticDislocations3d.normalFaultsPositionAtDeepFaultTip(xDeep,zPosition,yDeep,alongStrikeLength,depthExtent,angle,np.zeros_like(angle),slip)
        faultsDataFrame=faults.ReturnDataFrame()
        
        return faultsDataFrame
    
    def LoadFaults(self, faultsDataFrame):
        faults=elasticDislocations3d.LoadNormalFaultsAtMidPoint(faultsDataFrame)
        return faults
        
    
#%%

class longTermUpliftExtensionalFixedLengthOfDomain(longTermUpliftExtensional,longTermUpliftFixedLengthOfDomain):
    pass
    
#%%
def ComputeLongTermDisplacement(dipAngle,distanceLocked,distanceSemiLocked,cutoff,mesh=None,sampleSizeForOneRun=1e7,minRationOfMeanToStd=0.2):
    """ compute long term uplift for dipAngle,distanceLocked,distaneSemi locked and cutoff, can also use a mesh given from above.
    return xarrays object. Can also give sample size
    This version of the function runs until max std reaches value lower than cutoff value"""

    t1=time.time()
    if mesh is None:
        upperPlateMesh=GenerateMesh(dipAngle)
    else:
        upperPlateMesh=mesh
        
    smallWedge=miscClasses.wedgeGeo(theta=dipAngle,x=distanceLocked)
    bigWedge=miscClasses.wedgeGeo(theta=dipAngle,x=distanceLocked+distanceSemiLocked)

    upperPlateColumnb=ComputeNormalizedColumbStressChange(dipAngle,smallWedge.R,bigWedge.R-smallWedge.R,upperPlateMesh)
    interp_cubic_geom=ComputeInterpolation(upperPlateMesh,upperPlateColumnb)

    cutoffPolygon=GetCutoffReturnContourPolygon(upperPlateMesh,upperPlateColumnb,cutoff)
    
    reallyBigWedge=miscClasses.wedgeGeo(theta=30,h=np.max(cutoffPolygon[:,1]))
    maxMw=RetrunMwForDepthExtent(reallyBigWedge.R)
    
   
    #rupturePoly=SetRupturePolygon(cutoffPolygon)

    lengthOfDomain=5
    grid=elasticDislocations3d.grid(x=np.arange(0,253,4),y=np.linspace(0,lengthOfDomain,70))
    totalDisplacment=grid.GenerateEmptyXRGrid()
    ratioOfMeanToStd=1000
    
    c=0
    faultsList=[]
    while ratioOfMeanToStd>minRationOfMeanToStd and c<35:
        
        Mw=gutenbergRichter.RandomSampleForGR(b=1, minMw=5, maxMw=maxMw,sampleSize=sampleSizeForOneRun)
        xDeep,yDeep,angle,Mw,depthExtent=ReturnEarthQuckesPositionAndAngle(Mw,cutoffPolygon,cutoffPolygon,interp_cubic_geom)
        slip,alongStrikeLength=ComputeSlipAndLength(Mw,depthExtent)
        
        
        zPosition=GenerateRandomFloatBetweenMinAndMax(0,lengthOfDomain,len(angle))


        faults=elasticDislocations3d.thurstFaultsPositionAtDeepFaultTip(xDeep,zPosition,yDeep,alongStrikeLength,depthExtent,angle,np.zeros_like(angle),slip)
        faultsDataFrame=faults.ReturnDataFrame()
        faultsDataFrame['counter']=c
        faultsList.append(faultsDataFrame)
        
        #faults=elasticDislocations3d.thurstFaultsPositionAtShallowFaultTip(xDeep,zPosition,yDeep,faultLenthAlongStrike,depthExtent,angle,np.zeros_like(angle),slip)
        
        defomation=elasticDislocations3d.okadaFor3dFaults(grid,faults)
        displacment=defomation.ReturnDisplacemnt()
        
        totalDisplacment=displacment+totalDisplacment
        
        ratioOfMeanToStd=CheckStdToMeanRation(totalDisplacment.dz)
        print("adding more EQs because ratioOfMeanToStd is: "+str(np.floor(ratioOfMeanToStd*100)),flush='True')
        c=c+1
    
    
    
    
    allFaultsDataFrame=pd.concat(faultsList,ignore_index=True).reset_index(drop=True)
   
    
    print("Finsihed compuation in mins: " + str(np.floor((time.time()-t1)/60)),flush='True')
    
    return totalDisplacment,allFaultsDataFrame,ratioOfMeanToStd

def ComputeAndSaveLongTerm(path,dipAngle,distanceLocked,distanceSemiLocked,cutoff,counter,mesh=None,sampleSizeForOneRun=1e7,minRationOfMeanToStd=0.2):
    """ this function calls the compution functions ! and it saves everything"""
    print("dip angle "+str(dipAngle)+" length of Locked Fault "+str(distanceLocked)+ " length of width zone "+str(distanceSemiLocked),flush=True)
    try:
        displacment,EQs,std=ComputeLongTermDisplacement(dipAngle,distanceLocked,distanceSemiLocked,cutoff,mesh,sampleSizeForOneRun,minRationOfMeanToStd)
        
        filename="dipAngle_"+str(dipAngle)+"_distanceLocked_"+str(distanceLocked)+"_distanceToFullyCreeping_"+str(distanceLocked+distanceSemiLocked)+"_distanceSemiLocked_"+str(distanceSemiLocked)+"_cutoff_"+str(cutoff)+"_c_"+str(counter)+"_std_"+str(std)
        filenameDisp=path+"displacment-"+filename+".netcdf"
        filenameEQ=path+"EQs-"+filename+".feather"
        
        displacment.dz.to_netcdf(filenameDisp)
        EQs.to_feather(filenameEQ)
        
    except:
        print ("can't compute run number "+str(counter))
        
def SaveLongTerm(path,dipAngle,distanceLocked,distanceSemiLocked,cutoff,counter,EQs,displacment,std,alpha=0.2,mu=0.6,outerString=None):
    """ this function calls the compution functions ! and it saves everything"""
    print("saving dip angle "+str(dipAngle)+" length of Locked Fault "+str(distanceLocked)+ " length of width zone "+str(distanceSemiLocked),flush=True)
    
    distanceSemiLocked=np.round(distanceSemiLocked,2)
    distanceLocked=np.round(distanceLocked,2)
    dipAngle=np.round(dipAngle,2)
    distanceFullyCreeping=np.round(distanceLocked+distanceSemiLocked,2)
    std=np.round(std,2)
    alpha=np.round(alpha,2)
    mu=np.round(mu,2)
    cutoff=np.round(cutoff,2)
    
    try:
        filename="dipAngle_"+str(dipAngle)+"_distanceLocked_"+str(distanceLocked)+"_distanceToFullyCreeping_"+str(distanceFullyCreeping)+"_distanceSemiLocked_"+str(distanceSemiLocked)+"_cutoff_"+str(cutoff)+"_c_"+str(counter)+"_std_"+str(std)+"_mu_"+str(mu)+"_alpha_"+str(alpha)
        if outerString is not None:
            filename=str(outerString)+filename

        filenameDisp=path+"displacment-"+filename+".netcdf"
        filenameEQ=path+"EQs-"+filename+".feather"
        
        displacment.dz.to_netcdf(filenameDisp)
        EQs.to_feather(filenameEQ)
        
    except Exception as e:
        
        print ("can't save run number "+str(counter),flush=True)
        print ("The error is: "+str(e),flush=True) 
    


def GetPointsReturnMeshCutoffColumnbStress(xy,meshSize,faults,cutoff=0.05):

    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(xy,
        mesh_size=meshSize)
        mesh = geom.generate_mesh()

    meshForComputing=elasticDislocations3d.gridTriangular(mesh.points,mesh.cells[1].data)
    coulmbStressChange=ComputeNormalizedColumbStressChange(dipAngle='dontneed',lockingDistanceAlongFault='dontneed',
    LockingWidthAlongFault='dontneed',mesh=meshForComputing,alpha=0.2,mu=0.6,interfaceFaults=faults)
    coulmbStressChangeAbove5Per=np.where(coulmbStressChange>cutoff,coulmbStressChange,0)
    cutoffPolygon=GetCutoffReturnContourPolygon(meshForComputing,coulmbStressChangeAbove5Per,cutoff)

    return meshForComputing,coulmbStressChangeAbove5Per,coulmbStressChange,cutoffPolygon
        
        
    
    
