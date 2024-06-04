#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 20:15:10 2024

@author: bar
"""
#%%
import numpy as np
import matplotlib.pyplot as plt



#%%
class wedgeGeo:
    """ this is a simple class where I either give it h/R/x and it gives me everything else. I have to give it dip angle
    
    theta - dip angle
    x - distance along surface
    h - depth
    R - largest part of the trinagle
    
        x
    __________
    \ theta |
     \      |
      \     |
      R\    |
        \   | h
         \  |
          \ |
           \|
    
 """
    
    def __init__(self,theta,x=None,h=None,R=None):
        if theta<0 :
            raise ValueError(" Theta has to be >0 ")
            
        self.theta=theta
        self.c=np.cos(np.deg2rad(theta))
        self.s=np.sin(np.deg2rad(theta))
        self.t=np.tan(np.deg2rad(theta))
        
        if x is not None and x>0 and h is None and R is None:
            self.x=x
            self.ComputeWithX()
        
        elif R is not None and R>0 and h is None and x is None:
            self.R=R
            self.ComputeWithR()
            
        elif h is not None and h>0 and R is None and x is None:
            self.h=h
            self.ComputeWithH()
            
        else:
            raise ValueError("Please give me a positive dip angle (in deg) and only value that is not None for the sides of the wedge. I need it as a positive number please")
        
        
    def ComputeWithX(self):
        self.R=self.x/self.c
        self.h=self.x*self.t
        
    def ComputeWithR(self):
        self.h=self.s*self.R
        self.x=self.c*self.R
        
    def ComputeWithH(self):
        self.x=self.h/self.t
        self.R=self.h/self.s
        
    def Plot(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        ax.plot([0,self.x],[0,0])
        ax.plot([self.x,self.x],[0,self.h])
        ax.plot([0,self.x],[0,self.h])
        
        ax.text(self.x/2,0,"x="+str(self.x))
        ax.text(self.x,self.h/2,"h="+str(self.h))
        ax.text(self.x/2,self.h/2,"R="+str(self.R))
        ax.text(0,0,"theta="+str(self.theta))
        ax.set_aspect('equal')
        ax.invert_yaxis()
            
    
