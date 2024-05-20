#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 00:16:15 2021

@author: bar
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd


#%%
class gunterbergRichterFit:
    """ this package gets an array/list of Mw and then compute a,b for the Mw.
    It can also plot and cutoff EQs smaller than a ceritan min value"""
    
    def __init__(self,Mw,sigma=None):
       self.Mw,self.N=self.GetMWReturnN(Mw)
       self.Fit(sigma)
       self.sigma=sigma
        
           
    def GetMWReturnN(self,Mw):
        """ this function gets Mw  sort it from small to large and then number it from large to small """
        Mw=np.sort(Mw)
        N=np.flipud(np.arange(len(Mw)))
        return Mw,N
    
    
    def Fit(self,sigma=None):
        """ fit between Mw and N based on the GR function 
        please notice that if I use sigma then I fit the invertred form of GR so to add errors on y axis"""
        if sigma is None:
            popt, pcov = curve_fit(GR, self.Mw, self.N,absolute_sigma=True)
        else:
            popt, pcov = curve_fit(InvertedGR,  self.N[0:-1],self.Mw[0:-1],absolute_sigma=True,sigma=sigma[0:-1])
            
        self.a=popt[0]
        self.b=popt[1]
        self.pcov=pcov
        
    def Cutoff(self,Cutoff,down=True):
        """ remove Mw smaller then a cutoff value and then return a new class """
        #MwData=pd.DataFrame({'Mw':self.Mw})
        if down==True:
            #ind=MwData.loc[MwData.Mw>Cutoff].index
            mask=np.where(self.Mw>Cutoff,True,False)
            
        else:
            #ind=MwData=MwData.loc[MwData.Mw<Cutoff].index
            mask=np.where(self.Mw<Cutoff,True,False)
            
        #return gunterbergRichterFit(MwData.loc[ind,'Mw'].values,sigma=self.sigma[mask[)
        if self.sigma is None:
            return gunterbergRichterFit(Mw=self.Mw[mask])
        else:
            return gunterbergRichterFit(Mw=self.Mw[mask],sigma=self.sigma[mask])
    
    
    
    

        
    def PlotAll(self,ax=None):
        """ plot everything including a,b and fitted function and recored EQs"""
        if ax is None:
            fig,ax=plt.subplots()
            
        self.PlotData(ax=ax)
        self.PlotGRCurve(ax=ax)
        self.PlotText(ax=ax)

        ax.set_yscale('log')
        
    def PlotData(self,ax=None,**args):
        if ax is None:
            fig,ax=plt.subplots()
        
        ax.scatter(self.Mw,self.N,**args)
        
    def PlotGRCurve(self,ax=None,minMw=None,maxMw=None,**args):
        if ax is None:
            fig,ax=plt.subplots()
            
        if minMw is None:
            minMw=np.min(self.Mw)-0.2
        if maxMw is None:
            maxMw=np.max(self.Mw)+0.2
        x=np.linspace(minMw,maxMw ,1000)
        
        ax.plot(x,GR(x,self.a,self.b),**args)
        
        
    def PlotText(self,ax=None,**args):
        
        if ax is None:
            fig,ax=plt.subplots()
        
        ax.text(np.min(self.Mw),0,"a="+str(np.round(self.a,2)),**args)
        ax.text(np.max(self.Mw),np.max(self.N),"b="+str(np.round(self.b,2)),**args)
        

    def FitBasedOnBins(self,magnitudes, bin_width=0.1):
        """
        params magnitudes : numpy.array
        params bin_width : float
        
        returns a,b,bstd, n-values if above the earthquake count threshold
        else returns np.nans
        https://anu-rses-education.github.io/EMSC-2022/Notebooks/LAB-week8/LAB8-Gutenberg-Richter.html
        """
        length = magnitudes.shape[0]
        minimum = magnitudes.min()
        average = magnitudes.mean()
        b_value = (1 / (average - (minimum - (bin_width/2)))) * np.log10(np.exp(1))
        square_every_value = np.vectorize(lambda x: x**2)
        b_stddev = square_every_value((magnitudes - average).sum()) / (length * (length - 1))
        b_stddev = 2.3 * np.sqrt(b_stddev) * b_value**2
        a_value = np.log10(length) + b_value * minimum
        
        return a_value, b_value, b_stddev, length

#%% 
def GR(Mw,a,b):
    """ this is simply GR compution with a,b,and Mw """

    return 10**(a-b*Mw)

def InvertedGR(N,a,b):
    return (a-np.log10(N))/b

# def RandomSampleForGROld(a,b,minMw,maxMw,sampleSize=10000):
#     """ give me a,b for known GR and minMw and maxMw to generate sample of GR disurubation.
#     This is based on the idea that GR is already the CDF of the diusbtration. So I normalize it between Mw min and max (0->1 on both x and y axis)
#     and then invert it and generate sample between 0-1 and return radnom GR diusbtration
#     https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
#     This link was a great resource
#     I was stupid here and there is a bug - Bar 16th of May , 2023 - please use version below"""
    
#     dMw=maxMw-minMw
#     a=(a-minMw)/dMw
#     b=b*minMw/dMw
#     norm=GR(minMw,a,b)
#     offset=GR(maxMw,a,b)
#     N=np.random.rand(int(sampleSize))
#     return (a-np.log10(N*norm+offset))/b


# def RandomSampleForGRoldV2(a,b,minMw,maxMw,sampleSize=10000):
#     """ give me a,b for known GR and minMw and maxMw to generate sample of GR disurubation.
#     This is based on the idea that GR is already the CDF of the diusbtration. So I normalize it between Mw min and max (0->1 on both x and y axis)
#     and then invert it and generate sample between 0-1 and return radnom GR diusbtration
#     https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
#     This link was a great resource
#     this is the new and improved version
#     Bar 16th of May, 2023"""
    
#     #term0=b*np.log(10)
#     #term1=10**(a-b*(minMw+maxMw))
#     #termMin=(10**(b*minMw))
#     #termMax=(10**(b*maxMw))
#     #norm=term1*(termMin*(term0*(minMw-maxMw)-1)+termMax)/term0
#     #norm=((10**(a-b*(minMw+maxMw)))*((10**(b*minMw))*((b*np.log(10))*(minMw-maxMw)-1)+10**(b*maxMw)))/(b*np.log(10))
#     #norm=(10**(a-b*minMw)-10**(a-b*maxMw))/(b*np.log(10))
#     norm=(10**(a-b*minMw))/(b*np.log(10))
#     offset=GR(maxMw,a,b)
#     offset=0
#     N=np.random.rand(int(sampleSize))
#     (np.log10(N*norm+offset)-a)/-b
#     return (a-np.log10(N*norm+offset))/b

# def RandomSamplerGRWhichCutAtMax(b,minMw,maxMw,sampleSize=10000):
#     Mw=np.array([])
#     while len(Mw)<sampleSize:
#         Mw_i=RandomSamplerForGRNotTruncated(b,minMw,sampleSize)
#         mask=np.where(Mw_i>maxMw,False,True)
#         Mw=np.append(Mw,Mw_i[mask])
        
        
# #     return Mw[0:int(sampleSize)]
    

# def RandomSamplerForGRNotTruncated(b,minMw,maxMw=None,sampleSize=10000):
#     """ give me a,b for known GR and minMw and maxMw to generate sample of GR disurubation.
#     This is based on the idea that GR is already the CDF of the diusbtration. So I normalize it between Mw min and max (0->1 on both x and y axis)
#     and then invert it and generate sample between 0-1 and return radnom GR diusbtration
#     https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
#     This link was a great resource
#     this is the new and improved version
#     Bar 17th of May, 2023
#     I also used this website https://mxrap.com/js_docs/lib_Gutenberg-Richter.html#GRRelationship
#     The truncated CDF verion with both Mmin and Mmax
#     or this paper
#     https://www-sciencedirect-com.insu.bib.cnrs.fr/science/article/pii/S0029549310003237
#     The truncated CDF is also found here
#     https://core.ac.uk/download/pdf/162043867.pdf """
    
#     beta=b*np.log(10)
    
#     if  maxMw is None:
#         N=np.random.rand(int(sampleSize))
#     else:
#         Nmax=1-np.exp(-beta*(maxMw-minMw))
#         N=GenerateRandomFloatBetweenMinAndMax(0,Nmax,sampleSize)
    
#     #Mw=fs.GenerateRandomFloatBetweenMinAndMax(minMw, maxMw,N=sampleSize)
#     #CDF=(1-np.exp(beta*(maxMw-Mw)))/
#     Mw=minMw-(np.log(1-N))/beta
    
#     return Mw

def RandomSamplerForGR(b,minMw,maxMw=None,sampleSize=10000):
    """ give me b for known GR and minMw and maxMw to generate sample of GR disurubation.
    if maxMw is None doesnt cut it
    This is based on the idea that GR is already the CDF of the diusbtration. So I normalize it between Mw min and max (0->1 on both x and y axis)
    and then invert it and generate sample between 0-1 and return radnom GR diusbtration
    https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
    This link was a great resource
    this is the new and improved version
    Bar 17th of May, 2023
    I also used this website https://mxrap.com/js_docs/lib_Gutenberg-Richter.html#GRRelationship
    The truncated CDF verion with both Mmin and Mmax
    or this paper
    https://www-sciencedirect-com.insu.bib.cnrs.fr/science/article/pii/S0029549310003237
    The truncated CDF is also found here
    https://core.ac.uk/download/pdf/162043867.pdf 
    please notice a is not used"""
    
    beta=b*np.log(10)
    if maxMw is None:
        norm=1
    else:
        norm=(1-np.exp(-1*beta*(maxMw-minMw)))
    
    N=np.random.rand(int(sampleSize))

    Mw=minMw-(np.log(1-N*norm))/beta
    
    return Mw

def GenerateRandomFloatBetweenMinAndMax(minLim,maxLim,N=1):
        return minLim+(maxLim-minLim)*np.random.rand(N)
#%%
class gunterbergRichterFitOld:
    """ this package gets a vector totalEQsMw that contains Mw of EQs in speififc area within a time period (default 3 years). 
    It then bins (according to bins paramter) it and compute the a,b for the GR of the area
    you can find the retuernd GR disubruation in the object self.gunterbergRichterDistrubation """ 
    def __init__(self,totalEQsMw,timePeriod=3,bins=5):
        self.totalEQsMw=totalEQsMw
        self.bins=bins
        self.timePeriod=timePeriod
        self.gunterbergRichterDistrubation=self.FitAandB()
            
    def SumEvents(self,vector):
        sumevents=np.zeros_like(vector)
        
        for i in range(len(sumevents)):
            sumevents[i]=np.sum(vector[i:])
            
        return sumevents
    
        
    def FitAandB(self):    
        
        numOfEvents,self.binedMw=np.histogram(self.totalEQsMw,self.bins)   
        self.binedMw=self.binedMw[1:]         
        self.sumEvents=self.SumEvents(numOfEvents)        
        gunterbergRichterDist=self.FitGRDistrubtion(self.binedMw,self.sumEvents)
        
        return gunterbergRichterDist                            
            
    def GunterbergRichterDistrbution(self,Mw,a,b):
        return 10**(a-b*Mw)

    def FitGRDistrubtion(self,Mw,sumOfEvents):
        popt, pcov = curve_fit(self.GunterbergRichterDistrbution, Mw, sumOfEvents)
        a=popt[0]
        b=popt[1]
        
        return gunterbergRichterCompute(a=a,b=b)
      

    
    
    def PlotSumEvents(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
            
        ax.scatter(self.binedMw,self.sumEvents)
        
    
    def PlotHistogram(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)

   
        ax.hist(self.totalEQsMw,bins=self.bins)
#%%   
class gunterbergRichterCompute:
    def __init__(self,a,b,timePeriod=3,totalTime=200,Mw=None):
        self.a=a
        self.b=b
        self.timePeriod=timePeriod
        self.totalTime=totalTime
        if Mw is None:
            self.Mw=np.arange(1,9,0.5)
            
        
            
        self.ComputeEvents()
        
    def ChangeTimePeriodAndTotalTimeToComputeAgain(self,totalTime=200,timePeriod=3):
        self.timePeriod=timePeriod
        self.totalTime=totalTime
        numOfEvents=self.ComputeEvents()
        
        return numOfEvents
    
    
    def MwToMoment(self,Mw):
        return 10**(1.5*Mw+9.1)
    
    def MomentToMw(self,Moment):
        return (np.log10(Moment)-9.1)/1.5
        
    
    def ComputeEvents(self,Mw=None,totalTime=None):
        if Mw is None:
            Mw=self.Mw
            
        if totalTime is None:
            totalTime=self.totalTime
        
        self.sumOfEvents=self.GunterbergRichterDistrbution(self.a,self.b,Mw)*(totalTime/self.timePeriod)
        self.eventsPerMw=self.SubstractEvents(self.sumOfEvents)
        self.ComputeDataFrameForEvents(self.eventsPerMw)
        
        return self.eventsPerMw
    
    def ComputeDataFrameForEvents(self,numOfEvents):
        data={'numOfEvents':numOfEvents,'Mw':self.Mw}
        
        self.eventsDataFrame=pd.DataFrame(data)
        self.eventsDataFrame['Moment']=self.MwToMoment(self.eventsDataFrame['Mw'])
        
        
        
    def GunterbergRichterDistrbution(self,a,b,Mw):
        return 10**(a-b*Mw)
    
    def SubstractEvents(self,sumOfEvents):

        numOfEvents=np.zeros_like(sumOfEvents)
        for i in range(len(sumOfEvents)-1,1, -1):
            numOfEvents[i]=np.sum(sumOfEvents[i:-1])
         
        numOfEvents[-2]=sumOfEvents[-2]-sumOfEvents[-1]
        numOfEvents[-1]=0
        
        return numOfEvents


    def PlotSumOfEvents(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.plot(self.Mw,self.sumOfEvents)
        
    def PlotNumOfEvents(self,ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
            
        ax.plot(self.Mw,self.eventsPerMw)
        
    def MinCutOff(self,minCutOff):
        self.eventsDataFrame=self.eventsDataFrame.loc[self.eventsDataFrame.Mw>minCutOff]
        
    def MaxCutOff(self,MaxCutOff):
        self.eventsDataFrame=self.eventsDataFrame.loc[self.eventsDataFrame.Mw<MaxCutOff]
        
    def FloorAndAddSesimicEventsFromMissingMoment(self,returnLeftoverEnergyInEventSize=7):
        leftOverEnergy=np.modf(self.eventsDataFrame['numOfEvents'])[0]*self.MwToMoment(self.eventsDataFrame['Mw'])
        allLeftOverEnergy=np.sum(leftOverEnergy)
        
        
        MomentForMissingEQs=self.MwToMoment(returnLeftoverEnergyInEventSize)
        missingEQ=np.round(allLeftOverEnergy/MomentForMissingEQs)
        
        self.eventsDataFrame.loc[:,'numOfEvents']=np.modf(self.eventsDataFrame.loc[:,'numOfEvents'])[1]
        
        oldEventsForMissingEQ=self.eventsDataFrame.loc[self.eventsDataFrame.Mw==returnLeftoverEnergyInEventSize,'numOfEvents']
        self.eventsDataFrame.loc[self.eventsDataFrame.Mw==returnLeftoverEnergyInEventSize,'numOfEvents']=oldEventsForMissingEQ+missingEQ
        
