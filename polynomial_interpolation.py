# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:17:39 2016

@author: bwinters
"""
import numpy as np
from functools import lru_cache

def newton_poly(X,Y,pVal):
    """This program evaluates a polynomial where X is the set of knots or 
    nodes to be interpolated and Y is the set of values taken at each value of X
    pVal is the value the polynomial is to be interpolated at"""
    
    # Make sure the knots and the values have the same number of points    
    if len(X)!=len(Y):
        print('length mismatch error')
        return(False)
    # Make sure no knots are repeated
    Xs={i for i in X}
    if len(Xs)!=len(X):
        print('input values not distinct')
        return(False)
    # Find the coefficents of the Newton polynomial
    xRange=range(1,len(X))
    c=[]
    c.append(Y[0])
    for k in xRange:
        d=X[k]-X[k-1]
        u=c[k-1]
        nList=list(range(k-1))
        nList.reverse()
        for i in nList:
            u=u*(X[k]-X[i])+c[i]
            d=d*(X[k]-X[i])
        c.append((Y[k]-u)/d)
    # Evaluate the Newton polynomial at pVal    
    p=1
    SUM=0
    for i in range(len(c)):
        if i>0:
            p=p*(pVal-X[i-1])
        SUM=SUM+p*c[i]
    return(SUM)

def lagrange_poly(knots,values,pVal):
    """This evaluates a polynomial that interpolates using the data in
    knots and values."""
    # Make sure the knots and values match
    lk=len(knots)
    lv=len(values)
    if lk!=lv:
        print('knots and values of different lenghts')
        return(False)
    # Make sure no knots are repeated
    kSet={k for k in knots}
    lks=len(kSet)
    if lks!=lk:
        print('knots not distinct')
        return(False)
    # Evaluate Lagrange Polynomials    
    pairs=[]
    for i in range(lk):
        pairs=pairs+[(knots[i],values[i])]
        
    y=0
    SUM=0
    for p in pairs:
        x,y=p
        SUM=SUM+y*lagrange_cardinal(x,knots,pVal)       
    return(SUM)
    
    
def lagrange_cardinal(root,knots,value):
    if not(root in knots):
        print('Root must be in knots')
        return(False)
    myList=[r for r in knots if r!=root]
    p=1
    for r in myList:
        p=p*(value-r)/(root-r)      
    return(p)

def divided_difference_poly(knots, values, pVal):
    lk=len(knots)
    lv=len(values)
    if lk!=lv:
        print('knots and values of different lenghts')
        return(False)
    kSet={k for k in knots}
    lks=len(kSet)
    if lks!=lk:
        print('knots not distinct')
        return(False)    
    c=divided_difference_table(knots, values)
    vrange=range(lv)
    coeffs=[]
    for i in vrange:
        coeffs.append(c[0,i])
    y=coeffs[0]
    termRange=range(1,lk)
    p=1
    for i in termRange:
       p=p*(pVal-knots[i-1]) 
       y=y+coeffs[i]*p
    return(y)
    
def hermite_poly(filename,t):
    c,X=hermite_dd_table(filename)
    coeffs=list()
    N=c.shape[0]
    for i in range(N):
        coeffs.append(c[0,i])
    prod=1
    S=coeffs[0]
    for i in range(1,len(X)):
        prod=prod*(t-X[i-1])
        S=S+coeffs[i]*prod
    return S

    
 
    
def divided_difference_table(knots,values):
    """This produces a divided difference table for the data given
    in knots and values.""" 
    ######
    # Check if the knots and the values match and are consistant
    #####
    lk=len(knots)
    lv=len(values)
    if lk!=lv:
        print('knots and values of different lengths')
        return(False)
    kSet={k for k in knots}
    lks=len(kSet)
    if lks!=lk:
        print('knots not distinct')
        return(False)
    ##########
    #Intialize c
    ###########
    c=np.matrix(np.zeros((len(knots),len(knots))))
    for i in range(lv):
        c[i,0]=values[i]
    #############
    # Populate c
    ###########
    jrange=range(1,lv)
    for j in jrange:
        irange=range(lv-j)
        for i in irange:
            c[i,j]=(c[i+1,j-1]-c[i,j-1])/(knots[i+j]-knots[i])
    return(c)

def hermite_dd_table(filename):
    """filename is a  comma separated file where each line contains:
    x,f(x),f'(x),...  This data is loaded and converted to a matrix
    upon which the divided difference algorithm is applied."""
    
    # Load Data    
    mfile=open(filename,'r')    
    table=list()
    for line in mfile:
        row=line.split(',')
        table.append([float(e) for e in row])
    mfile.close()
    # Put data into a matrix c.
    X=list()
    Y=list()
    for row in table:
        X=X+(len(row)-1)*[row[0]]
        Y=Y+(len(row)-1)*[row[1]]
    none_mat=list()
    for i in range(len(X)):
        none_mat.append(len(X)*[None])
    c=np.matrix(none_mat)
    for i in range(len(Y)):
        c[i,0]=Y[i]
    level=0
    for row in table:
        for i in range(1,len(row)):
            c[level,i-1]=row[i]/factorial(i-1)
        level=level+len(row)-1
    
    ##########
    # Divided Difference Algorithm
    ##########
    N=len(X)
    for j in range(1,N):
        for i in range(N-j):
            if X[i+j]!=X[i]:
                c[i,j]=(c[i+1,j-1]-c[i,j-1])/(X[i+j]-X[i])
    return (c,X)

    
    
def chebyshev_nodes(a,b,n):
    """This returns n Chebyshev nodes in [a,b]."""
    mR=range(1,n-1)
    nodes=[]
    for k in mR:
        mNode=.5*(a+b)+.5*(b-a)*np.cos(np.pi*(2*k-1)/(2*n))
        #print(k,mNode)
        nodes.append(mNode)
    nodes.reverse()   
    return(nodes)

@lru_cache(maxsize=10000)   
def factorial(n):
    if n==0:
        return 1
    else:
        return (n*factorial(n-1))
        
def fill_out(r,n):
    s=r.copy()
    s=s+(n-len(s))*[0]
    return s