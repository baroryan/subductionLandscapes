#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 14:13:07 2022

@author: bar
"""
import numpy as np
import matplotlib.pyplot as plt
import elasticDislocations3d

#%%
class slipDistribution:
    
    """ this class is rather simple, desigen to get proproties of max and min slip disurbation and then to split based on certian
    computaion """
    
    def __init__(self,maxNormalizedSlip=1,minNormalizedSlip=0,numberOfFaultToSlip=10,maxPositionOfTheFault=1,minPositionOfTheFault=0):
        if minNormalizedSlip < 0 or maxNormalizedSlip > 1 :
            raise TypeError("position and slip rates are normalized and need to be betwwen 0-1")
        
        if minNormalizedSlip>=maxNormalizedSlip:
            raise TypeError("min has to be smaller than max")

        if minPositionOfTheFault < 0 or maxPositionOfTheFault > 1 :
            raise TypeError("position and slip rates are normalized and need to be betwwen 0-1")
        
        if minPositionOfTheFault>=maxPositionOfTheFault:
            raise TypeError("min has to be smaller than max")            
        
        
        self.maxNormalizedSlip=maxNormalizedSlip
        self.minNormalizedSlip=minNormalizedSlip
        self.numberOfFaultToSlip=numberOfFaultToSlip
        self.ComputeAllSlip(maxPositionOfTheFault, minPositionOfTheFault)
        #self.ComputeAllSlip()
        
    def ComputeAllSlip(self,maxPositionOfTheFault,minPositionOfTheFault):
        self.maxInd=self.numberOfFaultToSlip-int(np.floor((1-maxPositionOfTheFault)*self.numberOfFaultToSlip))
        self.minInd=int(np.floor(minPositionOfTheFault*self.numberOfFaultToSlip))
        self.slipDistribution=np.ones(self.numberOfFaultToSlip)
        
        self.slipDistribution[0:self.minInd]=np.ones(self.minInd)*self.minNormalizedSlip
        self.slipDistribution[self.maxInd:]=np.ones(self.numberOfFaultToSlip-self.maxInd)*self.maxNormalizedSlip
        
        numOfFaultsForVaryingFaults=self.numberOfFaultToSlip-(self.numberOfFaultToSlip-self.maxInd)-self.minInd
        slipDistributionWithVaryingCoupling=self.ComputeSlipDistribution(numOfFaultsForVaryingFaults)
        
        
        self.slipDistribution[self.minInd:self.maxInd]=slipDistributionWithVaryingCoupling
        
    def FlipSlipDistribution(self):
        self.slipDistribution=np.flipud(self.slipDistribution)
    
    #def ComputeAllSlip(self):
    #    self.ComputeSlipDistribution(self.numberOfFaultToSlip)
    
    # I finished doing this now I need to test it
    def ComputeSlipDistribution(self):
        pass
    
    def ReturnSlipDistriubation(self):
        return self.slipDistribution
#%%
class linearSlipDistribution(slipDistribution):
    """ this class inherts its proprites from slipDisubartion and uses linear disubartion """
    def ComputeSlipDistribution(self,numOfFaults):
        slipDistribution=np.linspace(self.minNormalizedSlip,self.maxNormalizedSlip,numOfFaults)
        return slipDistribution
        
class UniformSlip(slipDistribution):
    """ this class inherts its proprites from slipDisubartion and uses unifrom slip along the fault that is equal to maxNormalizedSlip """
    def ComputeSlipDistribution(self,numOfFaults):
        slipDistribution=np.ones(numOfFaults)*self.maxNormalizedSlip
        return slipDistribution
    
#class CrackTipDistribution
#%%
class splitDislocationAlongManySmallDisocationWithVaryingDipSlip:
    def __init__(self,fault,slipDistribution):
        
      self.CheckForFault(fault)  
      self.fault=fault
      self.slipDistribution=slipDistribution
      self.ComputerNewFaultsPosition()
        
    def ComputerNewFaultsPosition(self):
        
        numberOfNewFaults=len(self.slipDistribution.slipDistribution)
        LengthOfNewFaults=(self.fault.faultLengthAlongDip)/numberOfNewFaults
        x,z=self.fault.ReturnShallowAndDeepTipOfault()
        
        
        xNewFaults=np.linspace(np.min(x),np.max(x),numberOfNewFaults+1)[0:-1]
        zNewFaults=np.linspace(np.min(z),np.max(z),numberOfNewFaults+1)[0:-1]
        
        
        newSlip=(self.slipDistribution.slipDistribution)*self.fault.dipSlip
        
        yFault=np.ones_like(xNewFaults)*self.fault.yFault
        strikeAngle=np.ones_like(xNewFaults)*self.fault.strikeAngle
        dipAngle=np.ones_like(xNewFaults)*self.fault.dipAngle
       
        strikeSlip=np.ones_like(xNewFaults)*self.fault.strikeSlip
        tensileSlip=np.ones_like(xNewFaults)*self.fault.tensileSlip
        faultLengthAlongStrike=np.ones_like(xNewFaults)*self.fault.faultLengthAlongStrike
        LengthOfNewFaults=np.ones_like(xNewFaults)*LengthOfNewFaults
        
        
 
            
        self.newFaults=elasticDislocations3d.thurstFaultsPositionAtShallowFaultTip(xNewFaults,yFault,zNewFaults, 
                                                                    
                                                                    faultLengthAlongStrike,LengthOfNewFaults,
                                                                    dipAngle,strikeAngle,
                                                                    newSlip,strikeSlip,tensileSlip)
    def ReturnNewDislocationFaults(self):
        return self.newFaults
        
    def CheckForFault(self,fault):
        if len(fault.xFault) > 1:
            raise TypeError ("Please send one fault!")
            
            
    def PlotOriginalFault(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
            
        self.fault.PlotFaultsDepthView(ax)
    
    def PlotNewFaults(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        self.newFaults.PlotFaultsDepthView(ax)
        
#%%
def ConverListOfPointsToListOfdipAnglesAndLengths(xPoints,yPoints,lastFaultDipAngle=None,lastFaultLength=2000,slipDisrutbtionObbjects=None,dipSlips=None):
    """ this functions takes points of faults and convertes them to length and dipAngles and return subsequentDislocationsToBuildComplexGeometries object """
    
    if len(xPoints) != len(yPoints): #or len(xPoints)%2!=0:
        raise TypeError("xPoints and yPoints need to be at the same length and divded by 2")
        
    
        
    lengths=np.zeros(len(xPoints)-1)
    dipAngles=np.zeros_like(lengths)
    
    for i in range(len(xPoints)-1):
        dipAngles[i]=np.rad2deg(np.arctan((yPoints[i+1]-yPoints[i])/(xPoints[i+1]-xPoints[i])))
        lengths[i]=((yPoints[i+1]-yPoints[i])**2+(xPoints[i+1]-xPoints[i])**2)**0.5
        
        
    if lastFaultDipAngle is not None and lastFaultDipAngle > 0:
        dipAngles=np.append(dipAngles, lastFaultDipAngle)
        lengths=np.append(lengths, lastFaultLength)
        
    elif lastFaultLength is not None and lastFaultDipAngle<0:
        dipAngles=np.append(dipAngles, dipAngles[-1])
        lengths=np.append(lengths, lastFaultLength)
        
    
    return subsequentDislocationsToBuildComplexGeometries(dipAngles,lengths,yPoints[0],xPoints[0],slipDisrutbtionObbjects,dipSlips)
    
        

#%%
class subsequentDislocationsToBuildComplexGeometries:
    """ This class get the location of the first dislocation and then arrays of dip angles and length and slips and then builds subsuqent faults.
    This is intndeeded to do build complex geomtries for dip angles faults like ramp strctuers in the Himalyas.
    This class also gets a list of slip disurbtions objects so it split parts of the faults"""
    
    def __init__(self,dipAngles,lengthsOfDislocations,depthOfFirstDislocationShallowTip,distanceOfFirstDislocationShallowTip,slipDistributionObjects=None,
                 dipSlips=None,infiniteLength=2000):
        """ This function gets 
        dipAngles array
        lengthsOfDislocation array
        depthOfFirstLocation float
        distanceOfFirst dislocation float
        slipDisrutbtionObbjects list. This is a bit tricky but each object will operate on the corsspondeing fault in the array. 
        if you don't want to use just use None object [None,None,None] for list of 3 arrays for example. defult list of Nones
        dipSlips arrays of slip, default 1 for all dislocations
        infinteLength for faults . default 2000km"""
        
        self.dipAngles=dipAngles
        self.lengthsOfDislocations=lengthsOfDislocations
        self.infiniteLength=infiniteLength
        
        
        if dipSlips is None:
            dipSlips=np.ones_like(self.dipAngles)
        self.dipSlips=dipSlips
        
        if slipDistributionObjects is None:
            slipDistributionObjects=len(self.dipAngles)*[None]
        self.slipDistributionObjects=slipDistributionObjects    
        
        self.CheckIfArraysAreOfTheSameLength([self.dipAngles,self.lengthsOfDislocations,self.dipSlips,self.slipDistributionObjects])
        positionedDislocations=self.GenerateArrayOfDislocationsBasedOnLocation(depthOfFirstDislocationShallowTip,distanceOfFirstDislocationShallowTip)
        self.splittedDislocations=self.ApplySlipDistributionFaults(positionedDislocations,self.slipDistributionObjects)
        
        
    def ReturnSplittedDislocation(self):
        return self.splittedDislocations
    def CheckIfArraysAreOfTheSameLength(self,listOfObects):
        for i in range(len(listOfObects)-1):
            
            if len(listOfObects[i+1])!=len(listOfObects[i]):
                print(i)    
                raise TypeError(" Arrays/lists are not of the same length")
            
    def GenerateArrayOfDislocationsBasedOnLocation(self,depthOfFirstDislocationShallowTip,distanceOfFirstDislocationShallowTip):
        dislocations=[]
        for i in range(len(self.dipAngles)):
            if i==0:
                dislocation_i=self.GeneratetDislocation(depthOfDislocationAtShallowTip=depthOfFirstDislocationShallowTip,
                                                        distanceOfDislocationAtShallowTip=distanceOfFirstDislocationShallowTip,
                                                        length=self.lengthsOfDislocations[i],
                                                  dipSlip=self.dipSlips[i],dipAngle=self.dipAngles[i])
            else:
                distanceOfFault_i,depthOfFault_i=dislocation_i.ReturnDeepTipOfault()
                dislocation_i=self.GeneratetDislocation(depthOfDislocationAtShallowTip=depthOfFault_i,
                                                        distanceOfDislocationAtShallowTip=distanceOfFault_i,
                                                        length=self.lengthsOfDislocations[i],
                                                  dipSlip=self.dipSlips[i],dipAngle=self.dipAngles[i])       
                
            dislocations.append(dislocation_i)
            
        return dislocations
            
    def ApplySlipDistributionFaults(self,dislocations,slipDistributionObjects):
        """ please remember for None objects it does nothing for linearDisubrtion object in array self.slipDistributionObjects it splits the fault"""
        
        redistributedDislocations=[]
        
        for slipDistribution_i,dislocation_i in zip(slipDistributionObjects,dislocations):
            
            if isinstance(slipDistribution_i,slipDistribution) is True: ## if slipDistribution_i is slipDistribution, split the dislocation
            
                splitDislocation=splitDislocationAlongManySmallDisocationWithVaryingDipSlip(dislocation_i,slipDistribution_i)
                redistributedDislocations.append(splitDislocation.ReturnNewDislocationFaults())
                
            elif slipDistribution_i is None: ## if not do nothing
                redistributedDislocations.append(dislocation_i)
            else:
                raise TypeError("Slip disubrtion object has to be either slipDistribution or None if you prefer not to use it. This is neither")
                
         
        newDislocations=dislocation_i.MergeFaultsStrctureIntoOne(redistributedDislocations)    
        return newDislocations
       
        
    def GeneratetDislocation(self,depthOfDislocationAtShallowTip,distanceOfDislocationAtShallowTip,length,dipSlip,dipAngle):
        
        fault=elasticDislocations3d.thurstFaultsPositionAtShallowFaultTip(distanceOfDislocationAtShallowTip,0,depthOfDislocationAtShallowTip, 
                                                                    self.infiniteLength,length,dipAngle,0,dipSlip,strikeSlip=0,tensileSlip=0)
        
        return fault
        
        
