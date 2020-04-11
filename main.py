import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import sys

def Shape_Functions_1D(Xi,Degree,dFlag):
    if dFlag == 1:
        N = np.zeros((2,np.size(Xi)))
        N[0,:] = -0.5
        N[1,:] = 0.5
    else:
        N1 = 0.5*(1-Xi)
        N2 = 0.5*(1+Xi)
        N = np.array([N1,N2])        
    return N

def Quadrature(Degree):
    if Degree == 1:
        Xi = np.array([0])
        W = np.array([2])
    elif Degree == 2:
        Xi = np.array([-0.577350269189626,0.577350269189626])
        W = np.array([1, 1])
    return Xi,W              

def Mesh():
    dofs = n*2
    auxX = np.linspace(0,l,n+1)
    n1 = np.arange(0,n)
    n2 = np.arange(1,n+1)
    Topology = np.array([n1,n2])
    x = auxX[Topology.T.flatten()]
    Topology = np.arange(0,dofs).reshape((n,2)).T
    return Topology,x

def Initial_Conditions(dofs):
    u0 = np.zeros([dofs,1])
    return u0

def Boundary_Conditions(iStep):
    dofs = np.size(x)
    rDofs = 0
    u = np.zeros([dofs,1])
    u[rDofs] = np.sin((2*np.pi)*iStep/5);
    return u,rDofs

def Advection_Velocity(xyz):
    a = 1
    return a

def Runge_Kutta(u,iStep,incrTime):
    def Residual(u,iStep):
        def Operators(u0,ud):
            dofs = 2
            Xi,W = Quadrature(2)
            ng = np.size(Xi)
            N = Shape_Functions_1D(Xi,1,0)
            dN = Shape_Functions_1D(Xi,1,1)
            # Mass Matrix
            refM = np.zeros((dofs**2,ng))
            for ig in range(ng):
                iN = N[:,ig]
                refM[:,ig] = refM[:,ig] + (iN[:,None]*iN*W[ig]).flatten()
            M = np.zeros((dofs**2,n))
            for iElm in range(n):
                iNodes = Topology[:,iElm]
                iCoords = x[iNodes]
                Jacobian = dN.T@iCoords
                for ig in range(ng):
                    M[:,iElm] = M[:,iElm] + (Jacobian[ig]*refM[:,ig]).flatten() 
            # Convection Matrix
            Xi,W = Quadrature(1)
            ng = np.size(Xi)
            N = Shape_Functions_1D(Xi,1,0)
            dN = Shape_Functions_1D(Xi,1,1)        
            # Advection Velocity        
            refC = np.zeros((dofs**2,ng))
            for ig in range(ng):
                iN = N[:,ig]
                idN = dN[:,ig]
                refC[:,ig] = refC[:,ig] + (idN[:,None]*iN*W[ig]).flatten()
            C = np.zeros((dofs**2,n))    
            for iElm in range(n):
                iNodes = Topology[:,iElm]
                iCoords = x[iNodes]
                for ig in range(ng):
                    a = Advection_Velocity(N[:,ig]@iCoords)
                    C[:,iElm] = C[:,iElm] + (a*refC[:,ig]).flatten()    
            # Flow Matrix
            m = np.zeros((dofs,n))    
            for iElm in range(n):
                iNodes = Topology[:,iElm]
                iCoords = x[iNodes]
                a = Advection_Velocity(iCoords)
                normal = np.array([-1,1])
                am = 0.5*(a*normal-np.abs(a*normal))
                m[:,iElm] = m[:,iElm] + am
              
            return M,C,m
        dofs = np.size(u)
        # Boundary Conditions
        ud,rDofs = Boundary_Conditions(iStep)
        # Get Operators
        [M,C,m] = Operators(u,ud)
        R = np.zeros((dofs,1))
        for iElm in range(n):
            dofs = 2
            iDofs = Topology[:,iElm]
            iM = M[:,iElm].reshape((dofs,dofs)).T
            iC = C[:,iElm].reshape((dofs,dofs)).T
            im = m[:,iElm]
            iu = u[iDofs]
            iud = ud[iDofs]
            iuout = np.zeros((dofs,1))
            inc = np.array([-1,1])
            neighDofs = iDofs+inc
            index =  np.logical_and(im!=0,neighDofs>0,neighDofs<dofs*n)
            iuout[index] = u[neighDofs[index]]
            if rDofs in iDofs:
                index2 = rDofs==iDofs
                iuout[index2] = iud[index2]
            iFlow = im[:,None]*(iuout-iu)
            iR = -inv(iM)@(iC@iu+iFlow)
            R[iDofs] = iR
        return R
    R1 = Residual(u,iStep)
    R2 = Residual(u+R1*incrTime/2,iStep+incrTime/2)
    R3 = Residual(u+R2*incrTime/2,iStep+incrTime/2)
    R4 = Residual(u+R3*incrTime,iStep+incrTime)
    u = u + (R1+2*R2+2*R3+R4)*incrTime/6
    return u

if __name__ == '__main__':
    # User input
    fileName = sys.argv[0]
    l = float(sys.argv[1])
    n = int(sys.argv[2]) # number of elements
    T = float(sys.argv[3]) # span of time
    incrTime = float(sys.argv[4]) # steps in iterative temporal solver
    dofs = 2*n
    # Create Mesh
    Topology,x = Mesh()
    # Initial Conditions
    u = Initial_Conditions(dofs)
    # Initialize Results Matrix
    results = np.zeros((dofs,int(np.ceil(T/incrTime))))
    results[:,0] = u[:,0]
    # Temporal Integration
    steps = np.arange(0,T,incrTime) # One step less
    for index,iStep in enumerate(steps):
        u = Runge_Kutta(u,iStep,incrTime)
        results[:,index] = u[:,0]
    # Plot Results
    fig = plt.figure()  
    ax = fig.gca()
    plt.xlim(0,l)
    plt.ylim(-1,1)
    fig.show()
    for index,iStep in enumerate(steps):
        ax.clear()
        ax.plot(x,results[:,index])
        plt.pause(incrTime)
    plt.close(fig)
               
    