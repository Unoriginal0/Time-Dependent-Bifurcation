import matplotlib

import matplotlib.pyplot as plt
import math
from matplotlib.ticker import NullFormatter
import numpy
import sys
    
"""
Marko Kocic, 2020-08-16

This code accompanies the paper which analyses the mathematical model used when 
treating the topic of Arctic sea-ice in climate science. The main purpose of the 
code is to find the steady state solutions to the ordinary differential equation
which describes the Arctic sea-ice mathematically.
"""
    
    
def f(t):
    """
    Smooth Non-analytic function
    """
    if t <= 0:
        return 0
    else:
        return math.e**(-1/t)
        
def g(t):
    """
    Modified Smooth Non-analytic function used as a bump function
    """
    return f(t)/(f(t)+f(1-t))
    
def albedo(t):
    """
    Albedo using a heavily modified bump function
    """
    ice = 0.55
    water = 0.3
    meltpoint = 273.15
    variance = 10

    return (water-ice)*g((t-(meltpoint-variance))/(2*variance))+ice
    
def Secant(f,df,x0,x1,tol,iter):
    """
    Secant method for a function f, using a parameter for df.
    Tailored specifically for the model.
    Input: model, df parameter, x0 and x1 where the root might be, tolerance, amount of iterations
    Output: x where f(x)=0
    """
    xnm1 = x1
    xnm2 = x0    
    
    for i in range(0,iter):
        xn = xnm1 - f(xnm1,df)*(xnm1-xnm2)/(f(xnm1,df)-f(xnm2,df))
        xnm2 = xnm1
        xnm1 = xn
        
        if abs(f(xn,df))<tol:
            return xn
            
    print("Solution not found")
    return None

def extrema(ls):
    """
    Find extremas of a function using a list of values.
    Does not differentiate if minima or maxima.
    Input: Y of a function as a list
    Output: List of extremas
    """
    res=[]
    if ls[1]-ls[0]<0:
        delta = False
    else:
        delta = True
        
    for i in range(1,len(ls)-1):
        if not delta and ls[i+1]-ls[i]>0:
            res.append((i,i+1))
            delta = True
        elif delta and ls[i+1]-ls[i]<0:
            res.append((i,i+1))
            delta = False
    
    return res
    
def emptytuples(ls):
    """
    Creates a list of empty tuples for the bifurcation diagram.
    Input: List of values for bifurcation diagram
    Output: List of tuples
    """
    res = []
    for i in range(0,len(ls)+1):
        res.append([[],[]])
    return res
    
def rootguess(ls):
    """
    Given a list of Y-values, creates a list of possible locations for roots.
    Input: List of Y-values
    Output: List of intervals where roots lie
    """
    res = []
    for i in range(0,len(ls)-1):
        if ls[i]*ls[i+1] < 0:
            res.append((i,i+1))
    return res
    
def rootexact(rootls,df):
    """
    Given a list of intervals, guesses the exact root using the given df.
    Input: List of intervals, df as integer
    Output: List of roots from the given interval
    """
    res = []
    for el in rootls:
        res.append(Secant(model, df, el[0], el[1],10**(-9), 50))
    return res        
        
def model(T,df):
    """
    The climate model as a function of T and bifurcation parameter df.
    Input: Temperature T, parameter df.
    Output: The change of temperature 
    """
    S = 1360.8
    alpha = albedo(T)
    sigma = 5.670374419*10**(-8)
    epsilon = 0.612
    
    return (S*(1-alpha)/4-epsilon*sigma*T**4)+df
    
if __name__ == "__main__":
    # Between which intervals we want to plot the temperature change
    min = 200
    max = 370
        
    X = numpy.linspace(min,max,(max-min)*1000+1)
    Y = [model(T,0) for T in X]
    plt.plot(X,Y, 'b') 
    plt.grid(True)
    plt.ylabel("dT/dt")
    plt.xlabel("Temperature (K)")
    plt.show()
    
    # Extrema
    test = extrema(Y)
    test = [min+el[1]/1000 for el in test]
    
    
    # Bifurcation parameters
    deltaFmin = -30
    deltaFmax = 30
    deltaFarray = numpy.linspace(deltaFmin,deltaFmax,(deltaFmax-deltaFmin)*100+1)
    
    solutions = emptytuples(test)    
    
    # Grade A spaghetti, this part takes AGES to complete
    for df in deltaFarray:
        newY = [el+df for el in Y]
        roots = rootguess(newY)
        roots = [(min+el[0]/1000,min+el[1]/1000) for el in roots]
        exact = rootexact(roots,df)
        print(round((df-deltaFmin)/(deltaFmax-deltaFmin)*100, 1), end='\r')
        sys.stdout.flush()
        
        for y in exact:
            added = False
            for l in range(0,len(test)):
                if y < test[l]:
                    solutions[l][0].append(df)
                    solutions[l][1].append(y)
                    added = True
                    break
            
            if not added:
                solutions[len(test)][0].append(df)
                solutions[len(test)][1].append(y)
    
    # Writes the bifurcation diagram to a file, we do not want to do this every
    # time we want to plot the bifurcation diagram. Plotting is done in a
    # different program.
    f = open("notime.csv","w+")
    f.write("{0}\n".format(len(test))) # This part gives how many plots we have to make
    for el in solutions:
        for i in range(0,len(el[0])):
            x = el[0][i]
            y = el[1][i]
            f.write("{0},{1}\n".format(x,y))
        
        f.write("\n")
    
    f.close()
    
    print("Done!")