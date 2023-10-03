#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 12:12:32 2021

@author: bar
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches 
import pandas as pd
from matplotlib import cm
import imageio
import matplotlib as mpl
from matplotlib import ticker, cm
import scipy.io
from scipy import interpolate
from matplotlib.patches import Rectangle
from scipy import interpolate
from shapely import geometry 
import glob
from scipy.optimize import fsolve
import warnings
from scipy.interpolate import RegularGridInterpolator
import matplotlib.tri as mtri
#%%
class material:
    def __init__(self,density=2700,porePressure=0.7,friction=0.6,cohesion=0, possionRation=0,shearModuls=0,label=None):
        self.rho_s=density
        self.rho_e=density*(1-porePressure)
        self.porePressure=porePressure
        self.mu=friction     
        self.C=cohesion
        self.v=possionRation
        self.shearModuls=shearModuls
        self.youngsModuls=2*shearModuls*(1+self.shearModuls)
        self.label=label

#%%
class coupling:
    def __init__(self,couplingProfile,couplingPrfileDepth,lockingDepth,lockingDistance):
        self.couplingProfile=couplingProfile
        self.couplingPrfileDepth=couplingPrfileDepth
        self.L=lockingDistance
        self.h=lockingDepth
        self.ComputeLockingWidth()
        
    
        
    def ComputeDistanceFromTrenchForDepth(self,depth):
        p=np.polyfit([0,self.h],[0,self.L],deg=1)
        
        distanceFromTrench=np.polyval(p,depth)
        
        return distanceFromTrench
        
        
        
    def ReturnCreepingDistance(self):
        return self.L
        
    def ReturnTransationFraction(self):
        return self.frac
    
    def ReturnLockingDistance(self):
        return self.lockingTransation
    
    def ComputeLockingWidth(self):
        
        for i in range(len(self.couplingProfile)):
            if self.couplingProfile[i]<1:
                break
            
        
        self.lockingTransationDepth=self.couplingPrfileDepth[i-1]
        self.lockingTransation=self.ComputeDistanceFromTrenchForDepth(self.lockingTransationDepth)
        frac=self.lockingTransationDepth/self.h
        self.frac=frac
        
    def ReturnDistanceProfile(self):
        return self.ComputeDistanceFromTrenchForDepth(self.couplingPrfileDepth)

    
    def PlotAll(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        ax.scatter(self.couplingPrfileDepth,self.couplingProfile,color='r')
        ax.vlines(self.lockingTransationDepth,0,1)
        
    
        
            
#%%
class wedge:
    g=9.81
    def __init__(self,bulk,interface,dipAngle,surfaceAngle,L,trenchDepth=0,couplingProfile=None,couplingPrfileDepth=None,label=None): 
        self.dipAngle=np.deg2rad(dipAngle)
        self.surfaceAngle=surfaceAngle
        self.theta=np.deg2rad(dipAngle+surfaceAngle)
        self.alpha=np.deg2rad(surfaceAngle)
        self.L=L
        self.ComputeWedgeVertices()
        self.interface=interface
        self.bulk=bulk
        self.label=label
        self.couplingProfile=couplingProfile
        self.couplingPrfileDepth=couplingPrfileDepth
        self.trenchDepth=trenchDepth*1000
        if interface.porePressure < bulk.porePressure:
            raise ValueError(" interface's pore pressure is smaller than bulk pore pressure! ")
        self.interface.mu=self.interface.mu*(1-self.interface.porePressure)/(1-self.bulk.porePressure)
        
        self.GenerateShapelyObjectForWedge()
        self.couplingObject=coupling(self.couplingProfile,self.couplingPrfileDepth,self.h,self.L)


        
    def ComputeWedgeVertices(self):
        h=np.tan(self.dipAngle)*self.L
        L2=h/np.tan(-self.dipAngle+np.pi/2)
        self.L2=L2
        self.h=h
        
        
        # self.wedgeVertices=np.array([[0,0],
        #                           [self.L,-np.tan(self.alpha)*self.L],
        #                           [self.L,np.tan(self.dipAngle)*self.L]])
        
        self.wedgeVertices=np.array([[0,0],
                                  [self.L,h],
                                  [self.L+L2,0]])
        
        self.downgoingPlateVertices=np.array([[-50e3,0],
            [0,0],
                                  #[0,np.tan(self.dipAngle)*self.L],
                                  [self.L*10,np.tan(self.dipAngle)*self.L*10],[-50e3,np.tan(self.dipAngle)*self.L*10]])
        
        self.upperPlatePlateVertices=np.array([[self.L+L2,0],
                                  [self.L,h],
                                  [self.L+L2,h]])
        
    def InsideWedge(self,X,Y):
        xNew=X.copy()
        yNew=Y.copy()
        epslion=((Y[2,0]-Y[1,0])**2+(X[0,1]-X[0,2])**2)**0.5
        for i,x in enumerate(X[0,:]):
            yWedgeDown=x*np.tan(self.theta)
            for j,y in enumerate(Y[:,0]):
            
                
                if x<0 or x>self.L:
                    xNew[j,i]=np.nan
                    yNew[j,i]=np.nan
                    continue
                
               
                if y>yWedgeDown+epslion:
                    xNew[j,i]=np.nan
                    yNew[j,i]=np.nan
        
        return xNew,yNew
    
    def GenerateShapelyObjectForWedge(self):
        p=self.wedgeVertices/1000
        self.wedgeShapely=geometry.Polygon([[p[i][0], p[i][1]] for i in range(len(p))])

    
    def PlotWedgeBoundaries(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
        
        ax.add_patch(matplotlib.patches.Polygon(self.wedgeVertices/1000,edgecolor='magenta',facecolor='None'))
        
        return ax
    
    def PlotBlockDowngoingPlate(self,ax):
        ax.add_patch(matplotlib.patches.Polygon(self.downgoingPlateVertices/1000,edgecolor='None',facecolor='White',zorder=10))
       #a=1+1
        
    def PlotBlockUpperPlate(self,ax,zorder=11):
        ax.add_patch(matplotlib.patches.Polygon(self.upperPlatePlateVertices/1000,edgecolor='None',facecolor='White',zorder=zorder))
        
    def SetLimitsBasedOnWedge(self,ax):
        ax.set_xlim([0,(self.L2+self.L)/1000])
        ax.set_ylim([self.h/1000,0])
        
    def PlotOnlyWedge(self,ax):
        #ax.set_aspect(1)
        self.PlotBlockDowngoingPlate(ax)
        #self.PlotBlockUpperPlate(ax)
        #self.SetLimitsBasedOnWedge(ax)
        #self.PlotLowerBoundaryAsCoupling(ax)
        ax.invert_yaxis()
        #self.SetLimitsBasedOnWedge(ax)
        
    def PlotLowerBoundaryAsCoupling(self,x=None,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
                 
        x,y,couplingForPloting=self.ComputeCouplingForX()


        colorbar=ax.scatter(x/1000,y/1000,s=8,c=couplingForPloting,cmap='seismic',zorder=-10000)        
        
        return ax,colorbar
        
        
        
    def ComputeCouplingForX(self,x=None):
        if x is None:
            x=np.linspace(0,self.L,num=10000)
            
        p=np.polyfit([0,self.L],[0,self.h],deg=1)
        y=np.polyval(p,x)
        f = interpolate.interp1d(self.couplingPrfileDepth, self.couplingProfile)
        couplingForDisturbation=f(y)
        
        return x,y,couplingForDisturbation
        
        
            
        
            
        
    def PlotOnlyWedgeBoundaries(self,ax):
        ax.set_aspect(1)
        self.PlotBlockDowngoingPlate(ax)
        #self.PlotBlockUpperPlate(ax)
        self.SetLimitsBasedOnWedge(ax)
        #self.PlotLowerBoundaryAsCoupling(ax)
        
        ax.invert_yaxis()
        #self.SetLimitsBasedOnWedge(ax)
        

#%%
class stressesFunctions:
    
    def FindMaxOrienationForCoulombStress(self,mu=0.6):
        """ This function takes sigma_xx , sigma_yy and sigma_xy and compute the max
        coulmnb stress orenation. Meaning what is the orineantion of the fault that is most likely to break. 
        Here I'm assuming minus stresses is compression """

        angle=np.zeros_like(self.Sxy)
        #c=0
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        for i in range(len(self.Sxx[0,:])):
            for j in range(len(self.Sxx[:,0])):
                Sxx=self.Sxx[j,i]
                Syy=self.Syy[j,i]
                Sxy=self.Sxy[j,i]
                
                if np.isnan(Sxx) or np.isnan(Syy) or np.isnan(Sxy):
                    angle[j,i]=np.nan
                    
                else:
                
                    def Func(theta):
                        c=np.cos(theta)
                        s=np.sin(theta)
                        return c*Sxx+(c-mu*s)*(s*Sxy-c*Syy)+(mu*c+s)*(c*Sxy+s*Syy)-c*(-mu*c-s)*Sxy+s*(c-mu*s)*Sxy
                    
                    
                    root = fsolve(Func,np.arctan(mu), xtol=1e-06, maxfev=500)
                    #print("c="+str(c)+" i="+str(i)+" j="+str(j))
                    angle[j,i]=np.rad2deg(root[0])
                    #c=c+1
        angle=np.mod(angle,360)
        return angle
    
    
    def FindColumnbStressChangeAtAngle(self,angle,mu=0.6):
        
        Sxx=self.Sxx
        Syy=self.Syy
        Sxy=self.Sxy
        angle=np.deg2rad(angle)
        c=np.cos(angle)
        s=np.sin(angle)
           
        return Sxx*s-c*Sxy*(c-mu*s)+(Sxy*s-c*Syy)*(s+mu*c)
    
    def FindMaxColumnbStress(self,mu=0.6):
        """ this function compute the max Columnb sress based on the max orenation computed by function FindMaxOrienationForCoulombStress"""
        

        angle=self.FindMaxOrienationForCoulombStress(mu=mu)      
        return self.FindColumnbStressChangeAtAngle(angle,mu=mu)
    
    def ReturnInterpolation(self,arrayOfZ,xToInterpolate,yToInterpolate):
        """ This function gets an array of values that has to be identical in dim to self.X and self.Y
        and locations of points of x and y and return the interpolation at these points.
        Pleast note that x and y need to be in meters
        based on this link:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
            """

        #from scipy.interpolate import RegularGridInterpolator
        #X=CascadiaWedgesWithDiterich.listOfWedges[0].dstress.X[0,:]
        #Y=CascadiaWedgesWithDiterich.listOfWedges[0].dstress.Y[:,0]
        #interp = RegularGridInterpolator((Y,X), anglesForMaxColumnbStresses)

        #pts = np.array([[8e3, 222e3],
        #                [24e3, 5.2e3]])
        
        X=self.X[0,:]
        Y=self.Y[:,0]
        interp = RegularGridInterpolator((Y,X), arrayOfZ)
        pts=np.array([yToInterpolate,xToInterpolate])
        pts=np.transpose(pts)
        
        return interp(pts)
         
    def ComputeAllForStresses(self):
        self.ComputePrincipalStresses()
        self.ComputeShearMax()
        self.ComputePrincipalVectors()
    
    def ComputePrincipalStresses(self):
        avgStress=(self.Sxx+self.Syy)/2
        diveratricStress=(((self.Sxx-self.Syy)/2)**2+(self.Sxy)**2)**0.5
        
        self.S1=avgStress-diveratricStress
        self.S3=avgStress+diveratricStress
        
    def ComputeShearMax(self):
        self.shearMax=np.sqrt(0.25*((self.Sxx-self.Syy)**2)+self.Sxy**2)
    
    
        
    def ComputePrincipalVectors(self,Sxx=None,Syy=None,Sxz=None):
        """Compute Phi which is the angle between maximum stress Sigma 1 and x axis"""
        
        self.phi=0.5*np.arctan(2*self.Sxy/(self.Sxx-self.Syy))
        
    def ComputeRegionReachedYiedling(self,S1,S3,mu):
        print("b")
     
    def Plot2DwithColor(self,MatrixToPlot,ax=None,fig=None,normVector=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
            
            
        if normVector is None:    
            surf = ax.pcolor(self.X/1000, self.Y/1000, MatrixToPlot/1e6, cmap=cm.coolwarm,
                             linewidth=0, antialiased=False,shading='auto')
        else:
            try:
                norm = mpl.colors.Normalize(vmin=np.min(normVector), vmax=np.max(normVector))
            except TypeError:
                norm=normVector
            surf = ax.pcolor(self.X/1000, self.Y/1000, MatrixToPlot/1e6, cmap=cm.coolwarm,norm=norm,
                             linewidth=0, antialiased=False,shading='auto')
        if fig is not None:    
            fig.colorbar(surf)
        ax.invert_yaxis()
        ax.set_xlabel('Width [km]')
        ax.set_ylabel('Depth [km]')
        
        return ax
   
    
    def ComputeColumbFauilre(self,mu,cohesion):
        """ based on eq 8.40 in T&S. Will code my own equation next"""
        
        S1=-1*self.S1
        S3=-1*self.S3
        
        
        #regionReachedFauilre=np.zeros_like(S1)
        regionReachedFauilre=np.zeros_like(S1,dtype=bool)
        #regionReachedFauilreBoolean=np.zeros_like(S1,dtype=bool)
       
        muAmp=np.sqrt(1+mu**2)
        columnbMohrCretria=S1*(muAmp-mu)-S3*(muAmp+mu)-2*cohesion
        
        #regionReachedFauilre=np.where(columnbMohrCretria>0,np.nan,columnbMohrCretria)
        regionReachedFauilre=columnbMohrCretria
        self.regionReachedFauilre=regionReachedFauilre
        self.regionReachedFauilreBoolean=np.where(regionReachedFauilre>0,1,0)
        
        
    def ComputeColumbFauilreBar(self,mu,cohesion):
        phi=0.5*np.arctan(1/mu)
        
        Sigma_n=self.ComputeNormalStressOnPlane(phi)
        tau=self.ComputeShearStressOnPlane(phi)
        
        regionReachedFauilre=np.zeros_like(tau,dtype=bool)
        
        
        columnbMohrCretria=Sigma_n*mu+cohesion-tau
        regionReachedFauilre=(columnbMohrCretria >=0 )
        self.regionReachedFauilre=regionReachedFauilre
        
        
    def ComputeNormalStressOnPlane(self,phi):
        phi=np.deg2rad(phi)
        
        S1=-1*self.S1
        S3=-1*self.S3
        Sigma_n=(S1+S3)*0.5+0.5*(S1-S3)*np.cos(2*phi)
        
        return Sigma_n
    
    def ComputeShearStressOnPlane(self,phi):
        S1=-1*self.S1
        S3=-1*self.S3
        tau=-1*(S1-S3)*0.5*np.sin(2*phi)
        
        return tau
    
    def ComputeShearAndNormalStressesOnPlane(self,phi):
        sxx=-1*self.Sxx
        syy=-1*self.Syy
        sxy=-1*self.Sxy
        
        phi=np.deg2rad(phi)
        sigma_n=0.5*(sxx+syy)+0.5*(sxx-syy)*np.cos(2*phi)+sxy*np.sin(2*phi)
        tau=-0.5*(sxx-syy)*np.sin(2*phi)+sxy*np.cos(2*phi)
        
        return sigma_n,tau
        
        
    def PlotSyy(self,ax=None):
        
        ax=self.Plot2DwithColor(self.Syy,ax=ax)
        ax.set_title('Syy[MPa]')
        return ax
        
    def PlotSxy(self,ax=None):
        
        ax=self.Plot2DwithColor(self.Sxy,ax=ax)
        ax.set_title('Sxy[MPa]')
        return ax
    
    def PlotSxx(self,ax=None):
        
        ax=self.Plot2DwithColor(self.Sxx,ax=ax)
        ax.set_title('Sxx[MPa]')
        return ax
    
    def PlotS1(self,ax=None):
        ax=self.Plot2DwithColor(self.S1,ax=ax)
        ax.set_title('S1[MPa]')
        return ax

    def PlotS3(self,ax=None):
        ax=self.Plot2DwithColor(self.S3,ax=ax)
        ax.set_title('S3[MPa]')
        return ax
    
    def PlotCountoursWithColors(self,matrixToPlot,ax=None,countours=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
            
       
            
       # ax=self.Plot2DwithColor(matrixToPlot/1e6,ax=ax)
        
        if countours is not None:
            cs=ax.contourf(self.X/1000,self.Y/1000,matrixToPlot/1e6,levels=countours,linewidth=6,linestyle="--",locator=ticker.LogLocator())
            cs.clabel( colors='k',inline=True,fontsize=6)
        
        ax.invert_yaxis()
        return ax
    
    
    def PlotCountours(self,matrixToPlot,ax=None,countours=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
            ax.invert_yaxis()
            
       
            
       # ax=self.Plot2DwithColor(matrixToPlot/1e6,ax=ax)
        
        if countours is not None:
            cs=ax.contour(self.X/1000,self.Y/1000,matrixToPlot/1e6,levels=countours,linewidth=6,linestyle="--",locator=ticker.LogLocator())
            cs.clabel( colors='k',inline=True,fontsize=6)
        
        
        return ax
    
    def PlotArrowsForAngles(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
            
        ax.quiver(self.X/1000,self.Y/1000,np.cos(self.phi),np.sin(self.phi),scale_units='inches',scale=8)
        
        ax.invert_yaxis()
        ax.set_xlabel('Width [km]')
        ax.set_ylabel('Depth [km]')
        
    def PlotCountrf(self,arrayToPlot,ax=None,**args):
        if ax is None:
            fig, ax = plt.subplots()
            
        cs=ax.contourf(self.X/1000,self.Y/1000,arrayToPlot,**args)
        #fig.colorbar(cs)
        #cs.clabel( colors='k',fontsize=6)
        
        return ax,cs
    
    def PlotArrowsForStresses(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
            
        ax.quiver(self.X/1000,self.Y/1000,self.S1/1e6,self.S3/1e6)
        
        ax.invert_yaxis()
        ax.set_xlabel('Width [km]')
        ax.set_ylabel('Depth [km]')
        
    
    def PlotStressDifference(self,ax=None):
        ax=self.Plot2DwithColor(self.S1-self.S3,ax=ax)
        ax.set_title('S1-S3[MPa]')
        return ax
        
    def PlotStressOrinenation(self,ax=None):
        ax=self.Plot2DwithColor(np.rad2deg(self.phi)*1e6,ax=ax)
        ax.set_title('phi [degress]')
        return ax
    
    def PlotRegionReachingFaulire(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
            
        norm = mpl.colors.Normalize(vmin=3, vmax=8)
        
        plotMe=np.where(self.regionReachedFauilre>0,np.nan,self.regionReachedFauilre)
        plotMe=np.log10(-1*plotMe)
        X=self.X/1000
        Y=self.Y/1000
        surf=ax.pcolor(X,Y, plotMe, cmap='jet',
        linewidth=0, antialiased=False,shading='auto',norm=norm)
#        fig.colorbar(surf)
        ax.invert_yaxis()
        
        return ax
    
    
    def PlotDiamterOfReginReachingFauilre(self,ax=None,**args):
        if ax is None:
            fig, ax = plt.subplots()
            
       
        cs1=ax.contour(self.X/1000,self.Y/1000,self.regionReachedFauilre,levels=[0])
        #ax.clabel(cs1,**args)
        
#%%
class stressesFunctionsTriangles(stressesFunctions):
    def SaveTrianglurMesh(self,triMesh):
        self.triMesh=triMesh
    
    
    def Plot2DwithColor(self,MatrixToPlot,ax=None,fig=None,normVector=None,**args):
        if ax is None:
            fig, ax = plt.subplots()
            
        cb=ax.tricontourf(self.triMesh, MatrixToPlot,**args)
            
        ax.invert_yaxis()
        fig.colorbar(cb)
        ax.set_xlabel('Width [km]')
        ax.set_ylabel('Depth [km]')
        
        return ax

    def ReturnInterpolation(self,matrixToInterpolate):
        
        return mtri.CubicTriInterpolator(self.triMesh, matrixToInterpolate, kind='geom')
        
        
        

    
#%%
class stresses(stressesFunctions):
    def __init__(self,Sxx,Syy,X,Y,Sxy=None):

            self.Sxx=Sxx
            self.Syy=Syy
            self.X=X
            self.Y=Y
            if Sxy is None:
                self.Sxy=np.zeros_like(self.Sxx)
            else:
                self.Sxy=Sxy
                
#            self.CheckAllHaveTheSameShape([self.X,self.Y,self.Sxy,self.Sxx,self.Syy])
    
            
    def CheckAllHaveTheSameShape(self,listOfArrays):
        if len(listOfArrays[0].shape) < 2:
            raise TypeError("Some input arguments are not 2D")
            
        for i in range(len(listOfArrays)-1):
            if listOfArrays[i].shape != listOfArrays[i+1].shape:
                raise TypeError("Not all input arguments have the same shape")
        
    
    def __add__(self,other):
        Sxx=self.Sxx+other.Sxx
        Syy=self.Syy+other.Syy
        Sxy=self.Sxy+other.Sxy
        
        return stresses(Sxx=Sxx,Syy=Syy,X=self.X,Y=self.Y,Sxy=Sxy)
    
    def __mul__(self,multiplier):
        Sxx=self.Sxx*multiplier
        Syy=self.Syy*multiplier
        Sxy=self.Sxy*multiplier
        
        return stresses(Sxx=Sxx,Syy=Syy,X=self.X,Y=self.Y,Sxy=Sxy)
    
    
        

#%%
class stressesAdvanceWithSlip:
    def __init__(self,intialStress,dstress,ds,wedgeProprties,totalSlip,displacemenIinitial,slipRate):
        self.intialStress=intialStress
        self.dstress=dstress
        self.ds=ds
        self.wedgeProprties=wedgeProprties
        self.intialStress=intialStress
        self.totalSlip=totalSlip
        self.displacemenIinitial=displacemenIinitial
        self.slipRate=slipRate
        self.accumulationOfStressAndDisplacementPerYear=(self.slipRate/100)/self.ds     #ds is defined in meter and slip rate is defined in cm/yr
        self.AdvanceBaseOnSlip(totalSlip)
        
        
    def SaveOnlyLastTimeStep(self):
        
        self.stresses=[self.stresses[-1]]
        self.displacements=[self.displacements[-1]]
        
    def AdvanceBaseOnSlip(self,totalSlip):
        self.slip=self.CreateSlipArray(totalSlip)
        self.stresses=[]
        self.displacements=[]
        
        stress_i=self.intialStress

        self.fauilreForSlip=np.zeros_like(stress_i.X)
        
        for i in range(len(self.slip)):
            stress_i=stress_i+self.dstress

            stress_i.ComputePrincipalStresses()
            stress_i.ComputeColumbFauilre(self.wedgeProprties.bulk.mu,self.wedgeProprties.bulk.C) 
            self.fauilreForSlip=stress_i.regionReachedFauilreBoolean*self.ds+self.fauilreForSlip
            self.stresses.append(stress_i)
            self.displacements.append(self.displacemenIinitial*i)
            
            
        self.fauilreForSlip=np.where(self.fauilreForSlip>0,self.fauilreForSlip,np.nan)
                
        
    def CreateSlipArray(self,totalSlip):
        slip=np.arange(start=0,stop=totalSlip,step=self.ds)
        return slip
        
        
    def AdvanceBaseOnTime(self,dt,slipRateOfRegion,totalTime):
        
        timePerds=self.ds*(1/slipRateOfRegion)
        
        
        for i in range(int(totalTime/dt)):
            self.stresses=self.stresses+self.dstress*timePerds
            
        


    def ComputeFailureInTheBulk(self):
        self.stresses.ComputePrincipalStresses()
        self.stresses.ComputeColumbFauilre(self.wedgeProprties.bulk.mu,self.wedgeProprties.bulk.C)        


    def AdvanceBaseOnSlipIncermentAndShowWhatFails(self,totalSlip):
        
        self.stresses.ComputePrincipalStresses()
        
        self.fauilreForSlip=np.zeros_like(self.stresses.S1)
        for i in range(int(totalSlip/self.ds)):
            self.currectSlip=i*self.ds
            self.stresses=self.stresses+self.dstress     
            self.ComputeFailureInTheBulk()
            
            self.fauilreForSlip=self.stresses.regionReachedFauilreBoolean*self.ds+self.fauilreForSlip
            
        self.fauilreForSlip=np.where(self.fauilreForSlip>0,self.fauilreForSlip,np.nan)
   
    def PlotFailureForSlip(self,slip=None,ax=None,levels=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if slip is None:
            slipIndex=-1
            
        #self.stresses[slipIndex].Plot2DwithColor(self.fauilreForSlip*1e6,ax=ax)    
        if levels is None:
            ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForSlip,ax=ax)    
        else:
            ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForSlip,ax=ax,levels=levels)    
        #surf=ax.pcolor(self.stresses.X/1000,self.stresses.Y, self.fauilreForSlip, cmap='jet',linewidth=0, antialiased=False,shading='auto')#,norm=norm)
        self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
        self.wedgeProprties.PlotOnlyWedge(ax=ax)
        
        return ax,cs
    
    def PlotVerticalDisplacementForSlipForYear(self,ax=None,color='None'):
        if ax is None:
            fig, ax = plt.subplots()

        slipIndex=1
  
        ax.plot(self.displacements[slipIndex].X/1000,self.displacements[slipIndex].Syy*self.accumulationOfStressAndDisplacementPerYear,color=color)
        
        return ax
        
    
    
    def PlotFailureForAsPrecent(self,totalSlip=None,slip=None,ax=None,levels=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if slip is None:
            slipIndex=-1
            
        if totalSlip is None:
            totalSlip=self.totalSlip
        args={'vmin':0,'vmax':100}
        #self.stresses[slipIndex].Plot2DwithColor(self.fauilreForSlip*1e6,ax=ax)    
        if levels is None:
            ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForSlip*100/self.totalSlip,ax=ax,vmin=0,vmax=100)    
        else:
            ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForSlip*100/self.totalSlip,ax=ax,vmin=0,vmax=100,levels=levels)
        #surf=ax.pcolor(self.stresses.X/1000,self.stresses.Y, self.fauilreForSlip, cmap='jet',linewidth=0, antialiased=False,shading='auto')#,norm=norm)
        self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
        self.wedgeProprties.PlotOnlyWedge(ax=ax)
        #cs = clippedcolorbar(cs)
        
        return ax,cs
    
    def ReturnFailureContourAsPrecent(self,totalSlip=None,slipIndex=None,ax=None,levels=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if slipIndex is None:
            slipIndex=-1
            
        if totalSlip is None:
            totalSlip=self.totalSlip

        #self.stresses[slipIndex].Plot2DwithColor(self.fauilreForSlip*1e6,ax=ax)    
        if levels is None:
            ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForSlip*100/self.totalSlip,ax=ax,vmin=0,vmax=100)    
        else:
            ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForSlip*100/self.totalSlip,ax=ax,vmin=0,vmax=100,levels=levels)
            
        plt.close()
        
        if cs is None:
            raise ValueError ("Couldn't get contours! ")
        
        return cs
        

    def PlotFailureForTime(self,slipRate,slip=None,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if slip is None:
            slipIndex=-1
        ##units are in cm/yr
        slipRate=slipRate*(0.01)
        
        self.fauilreForTime=self.fauilreForSlip/slipRate  
        ax,cs=self.stresses[slipIndex].PlotCountrf(self.fauilreForTime,ax=ax)    
        #surf=ax.pcolor(self.stresses.X/1000,self.stresses.Y, self.fauilreForSlip, cmap='jet',linewidth=0, antialiased=False,shading='auto')#,norm=norm)
        self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
        self.wedgeProprties.PlotOnlyWedge(ax=ax)
        
        return ax,cs
        

        
    def AdvanceBaseOnSlipWithAnimation(self,totalSlip,path): # old function not in used I think
        figsForAnimations=[]
         
        for i in range(int(totalSlip/self.ds)):
            self.currectSlip=i*self.ds
            self.stresses=self.stresses+self.dstress     
            self.stresses.ComputePrincipalStresses()
            
            self.stresses.ComputeColumbFauilre(self.wedgeProprties.bulk.mu,self.wedgeProprties.bulk.C)
            figurePath=self.PlotFailure(self.currectSlip, i)
            
            figsForAnimations.append(imageio.imread(figurePath))
        
        
        imageio.mimwrite(path,figsForAnimations,fps=3)
            
#            ax=self.stresses.PlotRegionReachingFaulire()
            #self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
            
    def PlotFailure(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        self.stresses.PlotRegionReachingFaulire(ax=ax)
        self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
        
        return ax
    
    def PlotFailureForAnimation(self,currentSlip,i,path='/Users/bar/GoogleDrive/phd/research/Arthur/ArthursDropBox/MyPythonCode/animation/'):  # old function not in use I think
        fig, ax = plt.subplots()
        self.stresses.PlotRegionReachingFaulire(ax=ax)
        self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
        ax.set_aspect('equal',adjustable="datalim")
        currentSlip=np.round(currentSlip,2)
        figurePath=path+str(i)+".png"
        ax.set_title("slip:"+str(currentSlip)+" meters")
        
        fig.savefig(figurePath,dpi=300)
        plt.close('all')
        
        
        return figurePath
    
    def PlotFauilrePerimeter(self,ax=None,indexOfWedgeStress=-1,showOnlyWedge=True,**args):
        if ax is None:
            fig, ax = plt.subplots()
            
            
        self.stresses[indexOfWedgeStress].PlotDiamterOfReginReachingFauilre(ax=ax,**args)
        if showOnlyWedge == True:
            self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
            self.wedgeProprties.PlotOnlyWedgeBoundaries(ax=ax)
                
 
        return ax
    

#%%
class stressesAdvanceWithSlipWithDiterichFailure(stressesAdvanceWithSlip):
    def ReturnFailureContourAsPrecent(self,totalSlip=None,slipIndex=0,ax=None,levels=None):
        if ax is None:
            fig, ax = plt.subplots()
            

       
        mu=0.6
        angle=np.ones_like(self.dstress.X)*30
        maxF=self.dstress.FindColumnbStressChangeAtAngle(angle=angle,mu=mu)
        maxF=np.where(maxF>0,np.nan,maxF) # set nan values that are far from failure, meaning >0
        maxF=(maxF/np.nanmin(maxF))*100
        
        
        if levels is None:
            ax,cs=self.stresses[slipIndex].PlotCountrf(maxF,ax=ax,vmin=0,vmax=100)    
        else:
            ax,cs=self.stresses[slipIndex].PlotCountrf(maxF,ax=ax,vmin=0,vmax=100,levels=levels)
            
        plt.close()
        
        if cs is None:
            raise ValueError ("Couldn't get contours! ")
        
        return cs
    
    def PlotFailureForAsPrecent(self,totalSlip=None,slipIndex=0,ax=None,levels=None):
        if ax is None:
            fig, ax = plt.subplots()
            

       
        #mu=self.wedgeProprties.bulk.mu*(1-self.wedgeProprties.bulk.porePressure)
        mu=0.6
        angle=np.ones_like(self.dstress.X)*30
        maxF=self.dstress.FindColumnbStressChangeAtAngle(angle=angle,mu=mu)
        maxF=np.where(maxF>0,np.nan,maxF) # set nan values that are far from failure, meaning >0
        maxF=(maxF/np.nanmin(maxF))*100
        
        
        if levels is None:
            ax,cs=self.stresses[slipIndex].PlotCountrf(maxF,ax=ax,vmin=0,vmax=100)    
        else:
            ax,cs=self.stresses[slipIndex].PlotCountrf(maxF,ax=ax,vmin=0,vmax=100,levels=levels)
            
        self.wedgeProprties.PlotWedgeBoundaries(ax=ax)
        self.wedgeProprties.PlotOnlyWedge(ax=ax)
        
        if cs is None:
            raise ValueError ("Couldn't get contours! ")
        
        return ax,cs

    
#%%
class manyMatLabWedges:
    def __init__(self,path):
        self.filesToLoad=sorted(glob.glob(path))
        self.MatLabWedges=[]
        
        for file in self.filesToLoad:
            wedge=dataFromMatLabForWedge(file)
            self.MatLabWedges.append(wedge)
            
    def PlotAllCouplingProfiles(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
        for matlabWedge in self.MatLabWedges:
            matlabWedge.PlotCouplingWithDepth(ax=ax)
    
#%% stupid class as dataStructure

class dataFromMatLabForWedge:
    def __init__(self,matLabFilename):
        dataFromMatLab=scipy.io.loadmat(matLabFilename)
        try:
            self.X=dataFromMatLab['x']
        except KeyError:
            self.X=dataFromMatLab['X']
        try:
            self.Y=dataFromMatLab['z']
        except KeyError:
            self.Y=dataFromMatLab['Z']
            
        self.dipAngle=float(dataFromMatLab['PARAMS']['dip'][0])  
        self.surfaceAngle=float(dataFromMatLab['PARAMS']['surfaceAngle'][0])
        self.trenchDepth=float(dataFromMatLab['PARAMS']['trenchDepth'][0])
        self.lockingDepth=float(dataFromMatLab['PARAMS']['D'][0])  
        self.couplingProfile=dataFromMatLab['couplingProfile'][0]
        self.couplingProfileDepth=dataFromMatLab['couplingProfileDepth'][0]
        duz=dataFromMatLab['duz'][0,:]*1000
        self.verticalDisplacment=stresses(Sxx=np.zeros_like(duz),Syy=-1*duz,X=self.X[0,:],Y=self.Y[0,:],Sxy=np.zeros_like(duz))
        self.ds=0.1
        self.intersesmicStress=stresses(Sxx=dataFromMatLab['dsxx'],Syy=dataFromMatLab['dszz'],X=self.X,Y=self.Y,Sxy=dataFromMatLab['dsxz'])


    def PlotCouplingWithDepth(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        ax.scatter(self.couplingProfileDepth,self.couplingProfile)

        
#%%
class wedgeFauilreForVariousProprites:
    def __init__(self,fileProprtiesFromMatLab,interfacesProprites,bulksProprties,slipRate,totalSlip=10 ,label=None):
        
        self.allMatLabWedges=manyMatLabWedges(fileProprtiesFromMatLab)
        self.matlabProprties=self.allMatLabWedges.MatLabWedges
        self.interfacesProprites=interfacesProprites
        self.bulksProprties=bulksProprties
        self.label=label
        self.totalSlip=totalSlip
        self.slipRate=slipRate
        
        self.listOfWedges=self.ComputeFailureForDifferentWedges(self.matlabProprties,interfacesProprites,bulksProprties,totalSlip)
        #self.ComputeIntersesmicStressesBasedOnSlip(possibleIntersesmicSlip,self.listOfWedges)
        
    def ComputeIntersesmicStressesBasedOnSlip(self,slip,listOfWedges):
        for wedge_i in listOfWedges:
            wedge_i.AdvanceBaseOnSlip(slip)
            wedge_i.ComputeFailureInTheBulk()
            
            

    def ComputeFailureForDifferentWedges(self,manyMatlabProprties,interfacesProprites,bulksProprties,totalSlip):
        if isinstance(interfacesProprites, list) is False:
            interfacesProprites=[interfacesProprites]
        if isinstance(bulksProprties, list) is False:
            bulksProprties=[bulksProprties]  
        if isinstance(totalSlip, list) is False:
            totalSlip=[totalSlip]
        if isinstance(manyMatlabProprties, list) is False:
           manyMatlabProprties=[manyMatlabProprties] 
           
            
         
        listOfWedges=[]
        for matlabProprties in manyMatlabProprties:
            for interface in interfacesProprites:
                for bulk in bulksProprties:
                    for slip in totalSlip:
                        
                    
                        wedge_i=wedge(bulk=bulk, interface=interface, dipAngle=matlabProprties.dipAngle, surfaceAngle=matlabProprties.surfaceAngle,
                                      trenchDepth=matlabProprties.trenchDepth,
                                    L=matlabProprties.lockingDepth/np.tan(np.deg2rad(matlabProprties.dipAngle)),label=str(bulk.label)+str(interface.label)+str(slip),
                                    couplingPrfileDepth=matlabProprties.couplingProfileDepth,couplingProfile=matlabProprties.couplingProfile)
                                    
                        
                        intialStresses=AndersonianFailureStressesOnWedge(wedge_i,matlabProprties.X,matlabProprties.Y)
                        intialStresses.ComputePrincipalStresses()
            
                        
                        stresses=self.LoadWedges(intialStressesFromFun=intialStresses,intersesmicStress=matlabProprties.intersesmicStress,
                                                 ds=matlabProprties.ds,wedge_i=wedge_i,totalSlip=slip,displacemenIinitial=matlabProprties.verticalDisplacment,
                                                 slipRate=self.slipRate)
                   
                        #stresses=stressesAdvanceWithSlip(intialStress=intialStresses,dstress=matlabProprties.intersesmicStress,
                         #                                             ds=matlabProprties.ds,wedgeProprties=wedge_i,totalSlip=slip,
                          #                                            displacemenIinitial=matlabProprties.verticalDisplacment,slipRate=self.slipRate)
                        
                        listOfWedges.append(stresses)
                
        return listOfWedges
    
    def LoadWedges(self,intialStressesFromFun,intersesmicStress,ds,wedge_i,totalSlip,displacemenIinitial,slipRate):
        
        stresses=stressesAdvanceWithSlip(intialStress=intialStressesFromFun,dstress=intersesmicStress,
                                                                    ds=ds,wedgeProprties=wedge_i,totalSlip=totalSlip,
                                                                      displacemenIinitial=displacemenIinitial,slipRate=slipRate)
        
        return stresses
        
    
    
    def SaveOnlyLastTimeStep(self):
        
        for wedge_i in self.listOfWedges:
            wedge_i.SaveOnlyLastTimeStep()
            
    
    
    def PlotFailurePerTime(self,indexForWedge=None,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if indexForWedge is None:
            indexForWedge=np.array(0)
                

                
        self.listOfWedges[indexForWedge].PlotFailureForTime(slipRate=self.slipRate,ax=ax)
        #ax.set_title(self.listOfWedges[indexForWedge].wedgeProprties.label)
        
        return ax 
    
    def PlotFailurePerTimeForAllWedges(self):
        fig, ax = plt.subplots(len(self.listOfWedges),1)
        if len(self.listOfWedges)==1:
            ax=[ax]
        
        for i,wedge_i in enumerate(self.listOfWedges):
            self.PlotFailurePerTime(indexForWedge=i,ax=ax[i])
            
            
    def PlotFailureAsPrecent(self,indexForWedge=None,ax=None,levels=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if indexForWedge is None:
            indexForWedge=np.array(0)
        
        ax,cs=self.listOfWedges[indexForWedge].PlotFailureForAsPrecent(totalSlip=self.totalSlip,ax=ax,levels=levels)
#        fig.colorbar(cs,ticks=[0,100])
        
        
        return ax,cs
    
    def GetFailurePrecentContuors(self,indexForWedge=None,ax=None,levels=None):
        if indexForWedge is None:
            indexForWedge=np.array(0)
        
        ax,cs=self.listOfWedges[indexForWedge].PlotFailureForAsPrecent(totalSlip=self.totalSlip,ax=ax,levels=levels)
        
        #plt.close()
        
        return cs

    
    def PlotFauilrePerimeter(self,ax=None,listOfWedges=None,showOnlyWedge=True):
        if ax is None:
            fig, ax = plt.subplots()
            
        if listOfWedges is None:
            listOfWedges=self.listOfWedges
            
        for wedge_i in listOfWedges:
            wedge_i.PlotFauilrePerimeter(ax=ax,showOnlyWedge=showOnlyWedge)
            
            if showOnlyWedge is True:
                wedge_i.wedgeProprties.PlotOnlyWedge(ax=ax)
            
        return ax
    
    
#%%
class wedgeFauilreForVariousPropritesWithDiterich(wedgeFauilreForVariousProprites):
    def LoadWedges(self,intialStressesFromFun,intersesmicStress,ds,wedge_i,totalSlip,displacemenIinitial,slipRate):
        
        stresses=stressesAdvanceWithSlipWithDiterichFailure(intialStress=intialStressesFromFun,dstress=intersesmicStress,
                                                                    ds=ds,wedgeProprties=wedge_i,totalSlip=totalSlip,
                                                                      displacemenIinitial=displacemenIinitial,slipRate=slipRate)
        
        return stresses

        
#%% I don't think I am using this part anymore

class AndersonianFailureStressesOnWedge(stresses):
    """ This class is based on the simple  argument that everything is at yield on the interface plus andrhesion stress field """
    def __init__(self,wedgeProprities,X,Y):
        self.wedgeProprities=wedgeProprities
        self.X=X
        self.Y=Y
    
        self.Sxx,self.Syy,self.Sxy=self.ComputeStresses(wedgeProprities,wedgeProprities.bulk,wedgeProprities.interface,X,Y)
                
        
    def ComputeStresses(self,wedge,bulk,interface,X,Y):
        
        topoWeight=bulk.rho_e*wedge.g*X*np.tan(wedge.alpha)
        rockWeight=bulk.rho_e*wedge.g*Y
        
        waterWeight=self.ComputeWaterWeight(wedge)
        
        Syy=topoWeight+rockWeight+waterWeight
        sigmaCofficent,cohesionCofficent=self.ComputeRationBetweenSxxAndSyy(interface.mu)
        
        Sxx=cohesionCofficent*interface.C+sigmaCofficent*Syy
        Sxy=np.zeros_like(Syy)
        
        Sxx=-1*Sxx
        Syy=-1*Syy
            
        return Sxx,Syy,Sxy
    
    def ComputeWaterWeight(self,wedge):
        
        if wedge.trenchDepth==0:
            return 0
        
        return 0
        # did not account for water depth because this is supposedely in the pore pressure
        
        waterHeight=wedge.trenchDepth-self.X*np.tan(wedge.alpha)
        waterHeight=np.where(waterHeight>0,waterHeight,0)   ## if this is above surface then water weight is 0
        waterWeight=waterHeight*wedge.g*1000
        
        return waterWeight
        
        
        
        
    def ComputeRationBetweenSxxAndSyy(self,mu):
        phi=np.arctan(mu)
             
        csc=1/np.sin(phi)
        cot=1/np.tan(phi)
        
        sigmaCofficent=(1+csc)/(csc-1)  #this experssion goes with Sigma 3
        cohesionCofficent=2*(-cot+csc**2)/(csc-1) #this experssion is the one that is not realted to chonesion
        #return 2*np.sin(phi)/(1-np.sin(phi))
        return  sigmaCofficent,cohesionCofficent
#%%58
class strain:
    def ComputeAllForStrain(self):
        
        self.ComputeStrainBasedOnPlaneStrain()
        
    def ComputeStrainBasedOnPlaneStrain(self):
        v=self.wedgeProporties.v
        E=self.wedgeProporties.youngsModuls
        
        amp=(1+v)/E
        
        self.eps1=amp*(-1*self.S1*(1-v)+v*self.S3)
        self.eps3=amp*(-1*self.S3*(1-v)+v*self.S1)
        #the minus sign in the two places above is becuase I changed the geolgocal and enggereing stress convection so 
        # before comperssion is minus and now it is postive 
        
        
    def PlotEps1(self,ax=None):
        ax=self.Plot2DwithColor(self.eps1*1e6,ax=ax)
        ax.set_title('Eps 1[MPa]')
        return ax
        
    def PlotEps3(self,ax=None):
        ax=self.Plot2DwithColor(self.eps3*1e6,ax=ax)
        ax.set_title('Eps 1[MPa]')
        return ax        
    
            
            
#%%
class mohrCoulombYield:
    def __init__(self,stresses,effectiveFriction,cohesion):
        self.stresses=stresses
        self.effectiveFriction=effectiveFriction
        self.cohesion=cohesion
        
    def ComputeFailure(self):
        """ based on eq 8.40 in T&S. Will code my own equation next"""
        S1=self.stresses.S1
        S3=self.stresses.S3
        
        regionReachedFauilre=np.zeros_like(S1,dtype=bool)
        
        mu=self.effectiveFriction
        muAmp=np.sqrt(1+mu**2)
        columnbMohrCretria=S1*(muAmp-mu)-S3*(muAmp+mu)-2*self.cohesion
        
        regionReachedFauilre=(columnbMohrCretria >=0 )
        self.regionReachedFauilre=regionReachedFauilre
        
        return regionReachedFauilre
        
        
    def PlotRegionReachingFaulire(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            self.fig=fig
        
        X=self.stresses.X/1000
        Y=self.stresses.Y/1000
        surf=ax.pcolor(X,Y, self.regionReachedFauilre, cmap='Greys',
        linewidth=0, antialiased=False,shading='auto')
        fig.colorbar(surf)
        ax.invert_yaxis()            
        
        

#%%
class differtialStress:
    def __init__(self,dsxx,dsyy,perSlipUnit):
        self.dsxx=dsxx
        self.dsyy=dsyy
        self.perSlipUnit=perSlipUnit

        
    def ReturnRatioOfTotalSlipToSlipUnit(self,totalSlipToAdvanceFor):
        return (totalSlipToAdvanceFor/self.perSlipUnit)

    def AdvancePerSlip(self,totalSlipToAdvanceFor):
        
        ratio=self.ReturnRatioOfTotalSlipToSlipUnit(totalSlipToAdvanceFor)
        total_dsxx=self.dsxx*ratio  
        total_dsyy=self.dsyy*ratio
        
        return total_dsxx,total_dsyy    
#%%
class thermalForWedge:
    def __init__(self,thermalData):
        self.thermalData=thermalData
        
    def CorrectSoNoTopohraphy(self,thermalData=None):
        if thermalData is None:
            thermalData=self.thermalData
            
        surface=np.unique(thermalData['x[m]'])

        for surf in surface: ## go over on all possible x
            indexToUpdate=thermalData.loc[thermalData['x[m]']==surf].index  #for each x find the columna
            minValueForSurface=np.min(thermalData.loc[indexToUpdate]['y[m]'])    #find surface for that column
            thermalData.loc[indexToUpdate,['y[m]']]=thermalData.loc[indexToUpdate]['y[m]']-minValueForSurface
            
        return thermalData
    
    def PlotRawData(self,ax=None,s=2):
        if ax is None:
            fig, ax = plt.subplots()
        ax.scatter(self.thermalData['x[m]']/1000,self.thermalData['y[m]']/1000,c=self.thermalData['T[c]'],s=2)
        ax.invert_yaxis()
        
        return ax
        
    def ComputeInterpolation(self):
        self.f = interpolate.interp2d(self.thermalData['x[m]'], self.thermalData['y[m]'], self.thermalData['T[c]'], kind='cubic')
        
    def Interpolate(self,x,y):
        try:
            interpolateT=self.f(x,y)
        except AttributeError:
            self.ComputeInterpolation()
            interpolateT=self.Interpolate(x,y)
            
        return interpolateT
    
    def PlotInterpolationOnTopOfRealData(self,x,y,TformThermalModel=None,ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            
        if  TformThermalModel is None:
            TformThermalModel= self.f(x, y)
        
        self.PlotRawData(ax=ax)
        
        box=Rectangle([np.min(x)/1000,np.max(y)/1000],(np.max(x)-np.min(x))/1000,(-np.max(y)+np.min(y))/1000,linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(box)
        X,Y=np.meshgrid(x,y)
        ax.pcolor(X/1000,np.flipud(Y/1000),TformThermalModel)
        
                  
        return ax
    
    def ComputeInterpolationForWedge(self,wedge):
        X=wedge.listOfWedges[0].stresses.X
        Y=wedge.listOfWedges[0].stresses.Y
        
        x=X[0,:]
        y=Y[:,0]
        
        T=self.Interpolate(x,y)
        
        return X,Y,T
    
    def PlotInterpolationBasedOnWedge(self,wedge):
        X,Y,Tnew=self.ComputeInterpolationForWedge(wedge)
        x=X[0,:]
        y=Y[:,0]
        
        Tnew=self.CorrectForMinAndMaxValue(Tnew)
        self.PlotInterpolationOnTopOfRealData(x,y,TformThermalModel=Tnew)
        
    def CorrectForMinAndMaxValue(self,T,minValue=-10,maxValue=1440):
        Tnew=np.where( (T>0) & (T<1440),T,np.nan)
        
        return Tnew
        
        
        
        
#%%  
def clippedcolorbar(CS, **kwargs):
    from matplotlib.cm import ScalarMappable
    from numpy import arange, floor, ceil
    fig = CS.ax.get_figure()
    vmin = CS.get_clim()[0]
    vmax = CS.get_clim()[1]
    m = ScalarMappable(cmap=CS.get_cmap())
    m.set_array(CS.get_array())
    m.set_clim(CS.get_clim())
    step = CS.levels[1] - CS.levels[0]
    cliplower = CS.zmin<vmin
    clipupper = CS.zmax>vmax
    noextend = 'extend' in kwargs.keys() and kwargs['extend']=='neither'
    # set the colorbar boundaries
    boundaries = arange((floor(vmin/step)-1+1*(cliplower and noextend))*step, (ceil(vmax/step)+1-1*(clipupper and noextend))*step, step)
    kwargs['boundaries'] = boundaries
    # if the z-values are outside the colorbar range, add extend marker(s)
    # This behavior can be disabled by providing extend='neither' to the function call
    if not('extend' in kwargs.keys()) or kwargs['extend'] in ['min','max']:
        extend_min = cliplower or ( 'extend' in kwargs.keys() and kwargs['extend']=='min' )
        extend_max = clipupper or ( 'extend' in kwargs.keys() and kwargs['extend']=='max' )
        if extend_min and extend_max:
            kwargs['extend'] = 'both'
        elif extend_min:
            kwargs['extend'] = 'min'
        elif extend_max:
            kwargs['extend'] = 'max'
    return fig.colorbar(m, **kwargs)
            
        
        
        