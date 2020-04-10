import numpy as np
import matplotlib.pyplot as plt
import sys

def Shape_Functions_1D(Xi,Degree,dFlag):
    if dFlag == 1:
        dN1 = -0.5
        dN2 = 0.5
        N = np.array([dN1,dN2])
    else:
        N1 = -0.5*(1-Xi)
        N2 = 0.5*(1+Xi)
        N = np.array([N1,N2])        
    return N

def Quadrature(Degree):
    if Degree == 1:
        Xi = 0
        W = 2
    return Xi,W                

def Mesh(l,nElms):
    x = np.linspace(0,l,nElms+1)
    n1 = np.arange(0,nElms)
    n2 = np.arange(1,nElms+1)
    Topology = np.array([n1,n2])
    return Topology,x

def Initial_Conditions(nNodes):
    u0 = np.zeros([nNodes,1])
    return u0

def Boundary_Conditions():
    return u

def Operators(Topology,XYZ,u0,ud):
    return M,C,m

def Advection_Velocity():
    a = 1
    return a

def Runge_Kutta():
    # Boundary Conditions
    ud = Boundary_Conditions()
    # Advection Velocity
    a = Advection_Velocity()
    # Get Operators
    [M,C,m] = Operators(Topology,XYX,u0,ud,a)
    return u

if __name__ == '__main__':
    # User input
    fileName = sys.argv[0]
    l = sys.argv[1]
    n = sys.argv[2] # number of elements
    dofs = n+1
    time = sys.argv[3] # span of time
    incrTime = sys.argv[4] # steps in iterative temporal solver
    # Create Mesh
    Topology,x = Mesh(l,n)
    # Initial Conditions
    u = Initial_Conditions(dofs)
    # Initialize Results Matrix
    results = np.zeros([dofs,time/incrTime])
    results[:,0] = u
    # Temporal Integration
    steps = np.arange(0,time,incrTime) # One step less
    for index,iStep in enumerate(steps):
        u = Runge_Kutta(u,iStep)
        results[:,index] = u
        
    # Plot Results
        
        
    