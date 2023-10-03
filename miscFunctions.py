#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 23:10:15 2019

@author: bar
"""
import numpy as np
#import pandas
import pandas as pd
import glob
import random
##import miscClasses
import xarray as xr
import inpoly
import matplotlib.pyplot as plt
from numba import jit, njit
import numba
##import subductionSculpting as fs



#%%

def interactive_plot(fig, ax):
    """ this a function that chat GPT help me write. 
    it gets fig and axis and then get into intaective mode where left click save points
    and right click remove the point. space finish it and return a list of the coords"""
    coords = []
    ax.set_title('ineractive mode, left click add point, right click del the nearest one, space done')
    lines = ax.get_lines()
    assert len(lines) != 0, "No initial line found"

    def onclick(event):
        nonlocal coords
        nonlocal lines
        ix, iy = event.xdata, event.ydata
        if ix is None or iy is None:  # click outside the plot
            return

        if event.button == 3:  # right click
            # Find and remove the closest point
            if coords:
                closest_point = min(coords, key=lambda coord: (coord[0]-ix)**2 + (coord[1]-iy)**2)
                coords.remove(closest_point)
                for line in lines:
                    if line.get_xdata()[0] == closest_point[0] and line.get_ydata()[0] == closest_point[1]:
                        line.remove()
                        lines.remove(line)
                fig.canvas.draw()

        elif event.button == 1:  # left click
            # Add a new point
            coords.append((ix, iy))
            lines.append(ax.plot(ix, iy,'x',color='red')[0])
            fig.canvas.draw()

    def on_key(event):
        nonlocal coords
        if event.key == ' ':
            fig.canvas.mpl_disconnect(cid)
            print("Stopped picking points.")
            coords = np.array(coords)  # convert to numpy array
            print("Coordinates: ", coords)
            plt.close(fig)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('key_press_event', on_key)


    
    return coords
#%%

def interactiveClicksOnPlot(ax=None,fig=None):
    if ax is None:
        fig,ax=plt.subplots()
        
    coords = []

    ax.set_title("Interactive picker, press space to finish")
    
    def onclick(event):
        nonlocal coords
        ix, iy = event.xdata, event.ydata
        print('x = %d, y = %d'%(ix, iy))
        ax.plot(ix,iy,marker='x', color='red')
        fig.canvas.draw()
        # append the coordinates to the coords list
        coords.append((ix, iy))


        return coords

    def on_key(event):
        nonlocal coords
        if event.key == ' ':
            fig.canvas.mpl_disconnect(cid)
            print("Stopped picking points")
            

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.show()
    
    return coords
    

#%%
def GetRandomPointsInsidePolygon(polygon,sampleSize=1000):
    sampleSize=int(sampleSize)
    xmin=np.min(polygon[:,0]);xmax=np.max(polygon[:,0]);ymin=np.min(polygon[:,1]);ymax=np.max(polygon[:,1])
    xrand=fs.GenerateRandomFloatBetweenMinAndMax(xmin,xmax,sampleSize)
    yrand=fs.GenerateRandomFloatBetweenMinAndMax(ymin,ymax,sampleSize)
    insidePolygon,_=CheckIfXandYptsInsidePolygon(xrand, yrand, polygon)
    return xrand[insidePolygon],yrand[insidePolygon]
    
def GetASetNumberOfPointsInsidePolygon(polygon,N=10000):
    N=int(N)
    
    x,y=GetRandomPointsInsidePolygon(polygon,sampleSize=N)
    
    while len(x)<N:
        xMore,yMore=GetRandomPointsInsidePolygon(polygon,sampleSize=N)
        x=np.append(x, xMore)
        y=np.append(y, yMore)
        
    return x[0:N],y[0:N]

def GenerateUniformRandomCircle(R=1,xCenter=0,yCenter=0,N=int(1000)):
    theta=fs.GenerateRandomFloatBetweenMinAndMax(0,2*np.pi,N)
    r=np.sqrt(fs.GenerateRandomFloatBetweenMinAndMax(0,R,N))
    
    return xCenter+r*np.cos(theta),yCenter+r*np.sin(theta)
    

#%%

def GetSlipAreaReturnMw(area,slip,elasticModuls=30e9):
    return (np.log10(area*slip*elasticModuls)-9.049)/1.5

def GetSlipRutprenLenthDepthReturnMw(ruptureLength,ruptureDepth,slip,elasticModuls=30e9):
    return GetSlipAreaReturnMw(ruptureLength*ruptureDepth,slip,elasticModuls)

#%%

def CheckIfXandYptsInsidePolygon(x,y,polygon):
    """ this function basically uses CheckIfPointsAreInsidePolygon but get x and y instead of points """
    pts=np.zeros([len(x),2])

    pts[:,0]=list(x)
    pts[:,1]=list(y)
    pointsInsidePolygon,pointsOnTopOfBoundary=CheckIfPointsAreInsidePolygon(pts,polygon)
    #pointsInsidePolygon=is_inside_sm_parallel(pts,polygon)
    
    return pointsInsidePolygon,None

#%%



@jit(nopython=True)
def is_inside_sm(polygon, point):
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1  


@njit(parallel=True)
def is_inside_sm_parallel(points, polygon):
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean) 
    for i in numba.prange(ln):
        D[i] = is_inside_sm(polygon,points[i])
    return D  


def CheckIfPointsAreInsidePolygon(points,polygon):
    """ this is a simple wrapper for checking if points are inside a polygon using the inpoly package
    here I'm just checking that both numpy array of size [N,2] 
    you can find more methods here
    https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python"""
    
    points=np.asarray(points)
    polygon=np.asarray(polygon)
    
    if len(points.shape) != 2 or len(polygon.shape) != 2 or polygon.shape[1] !=2 or points.shape[1] != 2:
        raise TypeError("please make sure both to use 2D arrays! of size [N,2] ")
        
    if polygon.shape[0]<3:
        raise TypeError("polygon needs to include more points")
    
    if isinstance(points,np.ndarray) is False:
        raise TypeError("please make sure both are numpy arrays")
        
    if isinstance(polygon,np.ndarray) is False:
        raise TypeError("please make sure both are numpy arrays")
        
    try:
        pointsInsidePolygon,pointsOnTopOfBoundary=inpoly.inpoly2(points,polygon)
    except:
        pointsInsidePolygon=is_inside_sm_parallel(points,polygon)
        pointsOnTopOfBoundary=None
    
    return pointsInsidePolygon,pointsOnTopOfBoundary


#%%
def LoadDirWithXarrayFilesAndReturnAsDataFramce(path):
    files=glob.glob(path)
    
    data=miscClasses.getListWithUnderScoreAndValuesTranslateToTable(files)
    data=data.data
    dataFromArray=xr.open_dataarray(files[0])
    
    x=dataFromArray.x.values
    dz=np.zeros([len(data),len(x)])
    dz_std=np.zeros_like(dz)
    
    for i,file in enumerate(files):
        dataFromArray=xr.open_dataarray(file)
        dz[i,:]=dataFromArray.mean('y')
        dz_std[i,:]=dataFromArray.std('y')
        
        
    data['dz']=dz.tolist()
    data['dz_std']=dz_std.tolist()
    
    return data,x


#%%       
def ReturnArrayThatIsStoredAsAList(listOflists):
    data=np.zeros([len(listOflists),len(listOflists[0]) ])
    
    for i in range(len(data)):
        data[i,:]=listOflists[i]
                   
            
    return data     

def GenereateRandomIndexList(length,jump):
    index=[]
    for i in range(0,length,jump):
        index.append(random.randint(0,jump-1)+i)
        
    return index


#%%
def getObspyCatlogReturnDataFrameCatalog(obspyCatalog):
    times = []
    lats = []
    lons = []
    deps = []
    magnitudes = []
    magnitudestype = []
    focal_mechansim=[]
    dip=[]
    strike=[]
    rake=[]
    for event in obspyCatalog:
        if len(event.origins) != 0 and len(event.magnitudes) != 0:
            times.append(event.origins[0].time.datetime)
            lats.append(event.origins[0].latitude)
            lons.append(event.origins[0].longitude)
            deps.append(event.origins[0].depth/1e3)
            magnitudes.append(event.magnitudes[0].mag)
            magnitudestype.append(event.magnitudes[0].magnitude_type )
            strike_i=event.focal_mechanisms[0].nodal_planes.nodal_plane_1.strike
            dip_i=event.focal_mechanisms[0].nodal_planes.nodal_plane_1.dip
            rake_i=event.focal_mechanisms[0].nodal_planes.nodal_plane_1.rake
            focal_mechansim.append([strike_i,dip_i,rake_i])
            dip.append(dip_i)
            strike.append(strike_i)
            rake.append(rake_i)
            
    pdEvents = pd.DataFrame({'lat':lats,'long':lons,'depth':deps,
                       'mag':magnitudes,'type':magnitudestype,'focalMechansim':focal_mechansim,'time':times,'strike':strike,'dip':dip,'rake':rake})
    
    
    return pdEvents
                      


#%%
def FindInd(valueToFind,vector):
    vector=np.abs(vector-np.ones(len(vector))*valueToFind)
    Idx=vector.argmin()
    return Idx

#%%
def ReadTimeFromDissFlacFile(file):
        
        time = np.fromfile(file ,sep=" ")
        timeInYears=np.zeros(int(len(time)/2))
        c=0
        for o in range(len(time)):
            if o % 2 == 0:
                timeInYears[c]=time[o]/(1000*1000*3600*24*365.25)
                c=c+1
     
        return timeInYears    
    
#%% Function that might be useful
def ContactDataFramesIntoOneBasedOnFilePattern(path,filePattern):
    listOfDataFrames=[]
    files=glob.glob(path+filePattern)
    for file in files:
        listOfDataFrames.append(pd.read_csv(file))
        
    dataframe=pd.concat(listOfDataFrames,ignore_index=True)
    
    return dataframe
        
#%% read general file csv
def ReadCSVFile(filename,delimiter=None,header=None,skiprows=None):
    data=pd.read_csv(filename,delimiter=delimiter,header=header,skiprows=skiprows)
    return data

def ReturnGlobFiles(path,pattern):
    files=sorted(glob.glob(path+"/"+pattern))
    return files
    
def FitGPS(x,a,b,c,d):
        return a + b*x + c*np.sin(2*np.pi*x+d) 
    
def GetVpVsDensityReturnYoungusMoudlsPossionRatio(density,Vp,Vs):
    """ This function gets Velcotiy of S and P waves as well as denisty and returns Young Modlus and Poisson's ratio.
    It is based on this webpage https://subsurfwiki.org/wiki/P-wave_velocity """
    
    E=(density*(Vs**2)*(3*Vp**2-4*Vs**2))/(Vp**2-Vs**2) # Young's modlus
    v=(Vp**2-2*Vs**2)/(2*(Vp**2-Vs**2)) #Poisson's ratio
    
    
    return E,v


def GetVpVsDensityReturnShearModlusPossionRatio(density,Vp,Vs):
    """ This function gets Velcotiy of S and P waves as well as denisty and returns Shear Modulus and Poisson's ratio
    It is based on this webpage https://subsurfwiki.org/wiki/P-wave_velocity """
    
    E=density*(Vs**2)# Shear Modulus 
    v=(Vp**2-2*Vs**2)/(2*(Vp**2-Vs**2)) #Poisson's ratio
    
    
    return E,v


def GetDensityYoungusMoudlsPossionRatioReturnVpVs(density,E,v):
    """ just like functions above"""
    Vp=np.sqrt((E*(1-v))/(density*(1+v)*(1-2*v)))
    Vs=np.sqrt(E/(2*density*(1+v)))
    
    return Vp,Vs


def GetDensityShearModlusPoissonRatioReturnVpVs(density,mu,v):
    Vp=np.sqrt((2*mu*(1-v))/(density*(1-2*v)))
    Vs=np.sqrt(mu/density)
    
    return Vp,Vs

def GetShearModlusPossionRationReturnYoungModlus(mu,v):
    return 2*mu*(1+v)

def GetPossionRationShearModReturnLameFirstParamter(shearMod,v):
    return 2*shearMod*v/(1-2*v)

def GetPossionVpDensityReturnYoungModlusShearModlusLameFirstParamter(v,Vp,density):
    sqrt=np.sqrt((2*v/(1-2*v))+2)
    Vs=Vp/sqrt
    
    shearMod,v=GetVpVsDensityReturnShearModlusPossionRatio(density,Vp,Vs)
    youngMod,v=GetVpVsDensityReturnYoungusMoudlsPossionRatio(density,Vp,Vs)
    LameFirstParamter=GetPossionRationShearModReturnLameFirstParamter(shearMod,v)
    
    print ("in GPa: [shear mod,Yound mod, First Lame parameter]")
    return shearMod/1e9,youngMod/1e9,LameFirstParamter/1e9
    
    
#%%
def GenerateArrayOfRandomColors(lengthOfArray):
    
    
    color=["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])
       for j in range(lengthOfArray)]
    
    return color
    
        

#%%

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError ( "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError ( "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError ( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

#%% Convert for Mike
def ConvertFileWithBaseLineToDataFrame(file):
    f=open(file,'r')
    data=[]
    for line in f:
        string=line
      
        intialList=string.split(" ")
        listWithOutSpace=[i for i in intialList if i]
        day=listWithOutSpace[1].split(".")[1][0:-1]
        year=listWithOutSpace[1].split(".")[0]
        numericalDate=float(year)+float(day)/365
        listWithOutSpace.append(year)
        listWithOutSpace.append(day)
        listWithOutSpace.append(numericalDate)
        dataFrame=pd.DataFrame(data=listWithOutSpace)
        dataFrame=pd.DataFrame.transpose(dataFrame)
        data.append(dataFrame)
        
    finalData=pd.concat(data,ignore_index=True)
        
    return finalData
    
    
    
    
        
    pd.concat(data)
