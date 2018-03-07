# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:23:30 2017

@author: bwinters
"""

def graph(knots, values):
    """This graphs a spline with given knots and values"""
    import numpy as np
    import pylab as plt
    y=[]
    x=np.linspace(knots[0],knots[-1],1000)
    for t in x:
        y.append(spline_eval(t,knots,values))
    plt.plot(x,y)
    return None   

def spline_eval(x, knots,values):
    """This evaluates a cubic spline with given knots
    and corresponding values at the point"""
    mD=cubic_spline_coeffs(knots, values)
    i=get_interval_index(x,knots)
    A=mD['A'][i]
    B=mD['B'][i]
    C=mD['C'][i]
    D=mD['D'][i]    
    S=A*(knots[i+1]-x)**3+B*(x-knots[i])**3+C*(x-knots[i])+D*(knots[i+1]-x)
    return S
   
 

def get_z(knots,values):
    """This gets the z-values for a spline curve
    It is used in cubic_spline_coeffs"""
    n=len(knots)-1
    h=[]
    b=[]
    u=[0]
    v=[0]
    z=[0 for i in range(n+1)]
    for i in range(n):
        h.append(knots[i+1]-knots[i])
        b.append(6*(values[i+1]-values[i])/h[i])
    u.append(2*(h[0]+h[1]))
    v.append(b[1]-b[0])
    for i in range(2,n):
        u.append(2*(h[i]+h[i-1])-h[i-1]**2/u[i-1])
        v.append(b[i]-b[i-1]-h[i-1]*v[i-1]/u[i-1])
    z[n]=0
    for i in range(n-1,0,-1):
        z[i]=(v[i]-h[i]*z[i+1])/u[i]

    return(z,h,u,v)
    
def cubic_spline_coeffs(knots,values):
    """This gets the cubic spline coefficients
    for given knots and values
    It is used in spline_eval"""
    myData={}
    z,h,u,v=get_z(knots,values)
    myData['z']=z
    myData['h']=h
    myData['u']=u
    myData['v']=v
    n=len(knots)-1
    A=[]
    B=[]
    C=[]
    D=[]
    for i in range(n):
        A.append(z[i]/(6*h[i]))
        B.append(z[i+1]/(6*h[i]))
        C.append(values[i+1]/h[i]-z[i+1]*h[i]/6)
        D.append(values[i]/h[i]-z[i]*h[i]/6)
    myData['A']=A
    myData['B']=B
    myData['C']=C
    myData['D']=D
    return(myData)
    

def get_interval_index(x,knots):
    """This returns the index of the interval that
    contains x. It is used in spline_eval"""
    n=len(knots)-1
    if (x<knots[0]) or (x>knots[-1]):
        print('Not in Range')
        return None
    for i in range(n):
        #print(i,knots[i],x,knots[i+1])
        if (x>=knots[i]) and (x<=knots[i+1]):
            return i
    print('Should never get here')
    return None
