#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:41:53 2022

@author: bar
"""

from shapely import geometry 
from shapely import Point
import matplotlib.pyplot as plt
import numpy as np
from shapely.ops import cascaded_union
from shapely.ops import unary_union
import miscFunctions
import geopandas as gpd
from shapely.geometry import Polygon
from matplotlib import patches
#%%
class contourV2:
    def __init__(self,polys,levelValue,color=None):
        self.levelValue=levelValue
        self.color=color
        
        if isinstance(polys,list):
            polys=self.UnionPolys(polys)

          
        self.poly=polys
        self.GenerateGeopandaShape()
        
    def GenerateGeopandaShape(self):
        self.geopandaShape=gpd.GeoSeries(self.poly)
        
    def IntersectWithOtherShape(self,otherShape):
        self.poly=self.poly.intersection(otherShape)
        self.GenerateGeopandaShape()
        
    def UnionPolys(self,polys):
        return unary_union(polys)
    
    def PlotContour(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
                
        self.geopandaShape.plot(ax=ax,color=self.color,label=self.levelValue)
    
    
#%%
class getContourGenerateSheplyContours:
    def __init__(self,cs,intersectWithOtherShape=None):
        self.levelValue=cs.levels
        self.colorArray=cs.get_cmap()(np.linspace(0,1,len(self.levelValue)))
        self.CheckIfCSHoldsData(cs)
        self.GenerateLevelsOfContours(cs,intersectWithOtherShape)
        
    def CheckIfCSHoldsData(self,cs):
        for collection in cs.collections:
            path=collection.get_paths()
            
            if len(path) > 0:
                return 0
            elif len(path)==0:
                allPathsEmpty=True
            else:
                raise ValueError ("Serious porblem , not sure what to think of it")
                
        
        if allPathsEmpty is True:
            raise ValueError (" All Paths are empty")
                
            
    def GenerateLevelsOfContours(self,cs,intersectWithOtherShape=None):
        self.contours=[]
        colorArraytemp=[]
        levelValuetemp=[]
        self.shapeContourLevels=[]
        for i, collection in enumerate(cs.collections):
            #print ("collection "+str(i)+":")
            polys=[]    
            for j,path in enumerate(collection.get_paths()):
                if path.to_polygons():
             #       print("\t path "+ str(j)+":")
                    for npoly, polypoints in enumerate(path.to_polygons()):
                        x = polypoints[:, 0]
                        y = polypoints[:, 1]
              #          print("\t\t poly "+str(npoly)+":")
                            # be careful with the following---check to make sure
                            # your coordinate system expects lat first, lon second
                        poly_init = Polygon([coords for coords in \
                                                 zip( x,y)])
                        if poly_init.is_valid:
               #             print("valid!")
                            poly_clean = poly_init
                        else:
                #            print(" not valid!")
                            poly_clean = poly_init.buffer(0.)
                        if npoly == 0:
                            poly = poly_clean
                        else:
                            #poly = poly.difference(poly_clean) in the original code I found online this was the line of code. Howevere I changed it to this. At least for tri contouff
                            poly = poly.union(poly_clean) #new line
   
                        polys.append(poly)
                #else:
                 #   print("ER")
            contourLevel=contourV2(polys=polys,levelValue=self.levelValue[i],color=self.colorArray[i])
            
            if intersectWithOtherShape is not None:
                contourLevel.IntersectWithOtherShape(intersectWithOtherShape)
            
            if contourLevel.poly.is_empty is False:
                self.contours.append(contourLevel)
                ## only is the contour is not empty than I add the level and color, otherwise I don't need it
                levelValuetemp.append(self.levelValue[i])
                colorArraytemp.append(self.colorArray[i])
                self.shapeContourLevels.append(contourLevel.poly)  ## this is for poplatre Random EQ based on polygons package
                #print("here 10")
                
            #else:
             #   print("GG")
                
        self.levelValue=np.array(levelValuetemp)  ##regernareting the levelvalue and color array
        self.colorArray=np.array(colorArraytemp)
                
    def FindPointInWhatContourAndReturnPolygon(self,x,y):
        p=Point(x,y)
        
        for contour_i in self.contours:
            if contour_i.poly.contains(p):
                return contour_i.poly
            
        print ("Cant find a polygon containing the point",flush=True)
        return 0
        
                
    def PlotAllContours(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        for contour_i in self.contours:
            contour_i.PlotContour(ax=ax)
            
        return ax

    def UnionAllPolygons(self,listOfPolygons):
        union=listOfPolygons[0]
        for poly in listOfPolygons[1:]:
            union=poly.union(union)
            
        return union
    
    def ConvertPolyExteriorPoints(self,poly):
        if isinstance(poly, Polygon):
            return patches.Polygon(list(poly.exterior.coords))
        elif isinstance(poly, geometry.multipolygon.MultiPolygon):
            return patches.Polygon(list(poly.convex_hull.exterior.coords))
            #raise TypeError ("Can't deal with polygons that are Multipolygon! Please send me a regualr polygon of type shapely.geometry.polygon.Polygon")
        else:
           raise TypeError ("Not sure what are you trying to do to me but for me to give the permiter of a polygon I need type shapely.geometry.polygon.Polygon")
           
           


#%%
class oneContour:
    """ This is the lowest level obejct in the set of 3 classes. It basically desribe one perimter of a contour level
    it holds the color and levelvalue and the shaeply shape and verticies of the contours. It also have a bunch of plot routines """
    
    def __init__(self,verticies,levelValue,color=None):
        self.verticies=verticies
        self.levelValue=levelValue
        self.color=color
        self.shape=self.ConvertVertricesToShape(verticies)
        self.GenerateGeopandaShape()
    
    def GenerateGeopandaShape(self):
        self.geopandaShape=gpd.GeoSeries(self.shape)
        
    def ConvertVertricesToShape(self,verticies):
         
        poly_init = geometry.Polygon([[verticies[i][0], verticies[i][1]] for i in range(len(verticies))])
        if poly_init.is_valid:
            poly_clean = poly_init
        else:
            poly_clean = poly_init.buffer(0.)
        
        #if npoly == 0:
        #    poly = poly_clean
        #else:
        #    poly = poly.difference(poly_clean)
        
        # this if statment comes from this link:
        # https://gis.stackexchange.com/questions/99917/converting-matplotlib-contour-objects-to-shapely-objects            
            
        return poly_clean 
    
    def PlotShaeplyExterior(self,shape,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
        
        if isinstance(shape,geometry.multipolygon.MultiPolygon):

            for i,geom in enumerate(shape.geoms):
                if i==0:
                    ax.plot(*geom.exterior.xy,color=self.color,label=str(self.levelValue))
                else:
                    ax.plot(*geom.exterior.xy,color=self.color)
        elif isinstance(shape,geometry.polygon.Polygon):
            ax.plot(*shape.exterior.xy,color=self.color,label=str(self.levelValue))
            
    def PlotContour(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        self.geopandaShape.plot(ax=ax,color=self.color,label=self.levelValue)
            
    def PlotContourExterior(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        self.PlotShaeplyExterior(self.shape,ax=ax)
        
        return ax
    
    def IntersectWithOtherShape(self,otherShape):
        self.shape=self.shape.intersection(otherShape)
        self.GenerateGeopandaShape()
        
    
        
#%%            
class oneLevelContour:
    """ This class is the midlevel obeject for these 3 classes. It gets a contour level with a bunch of permiimters and generate objects for
    each permtier giving them level value and color. From here you can call the lower level plot functions"""
    
    def __init__(self,contourLevels,levelValue,color='None'):
        self.contourLevels=contourLevels
        self.levelValue=levelValue
        self.color=color
        self.GenerateShaply(contourLevels)
        self.OneShapeToRuleThemAll()
        #self.CleanHoles()
    
    
    def CleanHoles(self):
        for i,contour_i in enumerate(self.contours):
            shape_i=contour_i.shape
            if i == 0:
                contourMain=shape_i
            else:
                contourMain=contourMain.difference(shape_i)
        if len(self.contours)>0:
            self.contours=[self.contours[0]]
            self.contours[0].shape=contourMain
        
            
            
    def OneShapeToRuleThemAll(self):
        tempListOfShapes=[]
        for contour in self.contours:
            tempListOfShapes.append(contour.shape)
            
        #self.oneShapeForLevel=cascaded_union(tempListOfShapes)
        self.oneShapeForLevel=unary_union(tempListOfShapes)
        self.geopandaShape=gpd.GeoSeries(self.oneShapeForLevel)               
    
    def GenerateShaply(self,contourLevels):
        self.contours=[]
        for i in range(len(contourLevels)):
            p = contourLevels[i]
            if len(p) >3 : ## if don't have 3 points keep going
                contour=oneContour(p,self.levelValue,self.color)
                self.contours.append(contour)
            
        return self.contours
    
    def PlotContourLevel(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        self.geopandaShape.plot(ax=ax,color=self.color,label=self.levelValue)
        
        return ax
    
    def PlotContourLevelExterior(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        for contour in self.contours:
            contour.PlotContourExterior(ax=ax)
            
            
        return ax
            
    def PlotOneContourLevel(self,ax=None):
        self.contours.PlotShaeply(self.oneShapeForLevel)
        
    def IntersectWithOtherShape(self,otherShape):
        for contour in self.contours:
            contour.IntersectWithOtherShape(otherShape)
            
        self.OneShapeToRuleThemAll()
           
#%%
class allLevelsContours:
    """ This is the highest level object. It gets countour array and then using contour.allsegs it genreates the two lower level obects"""

    def __init__(self,cs,colorArray=None):
        self.cs=cs
        self.levelValue=cs.levels
        if colorArray is None:
            #self.colorArray=miscFunctions.GenerateArrayOfRandomColors(len(self.levelValue))
            self.colorArray=cs.get_cmap()(np.linspace(0,1,len(self.levelValue)))
        self.GenerateLevelsOfContours(cs)
            
    def GenerateLevelsOfContours(self,cs=None):
        if cs is None:
            cs=self.cs
            
        self.contoursLevels=[]
        self.shapeContourLevels=[]
        for i in range(len(cs.allsegs[:])):
            contoursForValue=oneLevelContour(cs.allsegs[i][:],levelValue=self.levelValue[i],color=self.colorArray[i])
            self.contoursLevels.append(contoursForValue)
            self.shapeContourLevels.append(contoursForValue.oneShapeForLevel)
            
    def PlotAllContoursExterior(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        for contourLevel in self.contoursLevels:
            contourLevel.PlotContourLevelExterior(ax=ax)
            
        return ax
    
    def PlotAllContours(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots()
            
        for contourLevel in self.contoursLevels:
            contourLevel.PlotContourLevel(ax=ax)
            
        return ax
    
    
    def IntersectWithOtherShape(self,shapeToIntersectWith):
        for contourLevel in self.contoursLevels:
            contourLevel.IntersectWithOtherShape(shapeToIntersectWith)
            
        
        
            
            
            
        
        
            
    
            
            
        
            
