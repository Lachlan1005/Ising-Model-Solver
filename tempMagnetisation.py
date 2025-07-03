import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import imageio
import random as rd
import os
import twoDIsing as tDI


def MagnetisationMetropolisFlips(N:int,J:float, T:float, maxIters:int=88000, verbosity:bool=False, lattice:list[list[int]] | None =None):
    """
    Return the magnetisation of the lattice in its thermal equillibrium state under temperature T and coupling strength J. 
    Solve using the Monte Carlo method, implementing the Metropolis model of interaction* 
    
    *: A spin flip is guaranteed to happen when it causes the overall energy to be lower. It can also happen if it raises the energy with a probabilty of P=exp(dE/T) where dE is the 
    energy raised by the flip

    Note: kB=1, express T in units of the coupling strength J
    """
    energies=[]
    iterationList=[]
    if lattice==None:
        lattice=tDI.gridGen(N)
    energyStable=0 
    iters=0
    while energyStable<=160 and iters<=maxIters:
        i=rd.randint(0,N-1)
        j=rd.randint(0,N-1)
        if verbosity:
            print("Trying to flip ", (j,i))
        lattice[i][j]=-lattice[i][j]
        newEnergy= tDI.hamiltonian(N,J,lattice)
        energies.append(newEnergy)
        iterationList.append(iters)
        f = (lattice[(i+1)%N][j] +lattice[(i-1)%N][j] +lattice[i][(j+1)%N] +lattice[i][(j-1)%N])
        dE=-2*J*lattice[i][j]*f
        if verbosity:
            print("Energy change if successsful:", dE)
        if dE>0:
            b=1/T
            transProb=np.exp(-b*dE)
            r=rd.random()
            if r>transProb:
                lattice[i][j]=-lattice[i][j] #flip back the lattice to original
                dE=0
                energyStable+=1
                if verbosity:
                    print("Flip failed")
            else:
                energyStable=0
                if verbosity:
                    print("Flip successful")
        elif dE<=0:
            energyStable=0
            if verbosity:
                print("Flip successful")
        if verbosity:
            for line in lattice:
                print(line)
      #  print(30*"\n","Solver running. Do not close or interrupt program.\n Iteration ", iters, "out of a maximum possibility of ", maxIters,"\n Equilibrium not reached yet. Current energy=", newEnergy, "Current Temperature=", T)
        if verbosity:
            print("Energy change from flip: ", dE)
            print("Energy: ", newEnergy)
            print("No Change Counter: ", energyStable)
   #     print(30*"=")
        iters+=1
    return abs(tDI.magnetisation(N, lattice))

def magnetisationTemperature(N:int,J:float,Tmin:float,Tmax:float,dT:float):
    """
    Plot magnetisation vs temperature for an NxN lattice
    Very memory intensive. Do not exceed N=5. 

    Comment out any calls to functions outside of the functions in twoDIsing.py before running this function
    """
    T=Tmin
    i=0
    temps=[]
    mags=[]
    while T<=Tmax:
        temps.append(T)
        mags.append(MagnetisationMetropolisFlips(N,J,T))
        print(30*"\n", "Solver running...")
        print("Current Temperature=", T, "|| Maximum Temperature=", Tmax)
        T+=dT
        i+=1
    print(30*"\n")
    print("Solving complete. See output plot for results.")
    plt.plot(temps,mags,color="black")
    plt.title("Magnetisation vs Temperature, N="+str(N))
    plt.xlabel("Temperature")
    plt.ylabel("Magnetisation")
    plt.show()
    return temps, mags

def terminalMagnetisationTemp():
    print(100*"=")
    print("Welcome to the Magnetisation-Temperature Plotter")
    N=int(input("Enter the size (N) of the NxN lattice: "))
    J=float(input("Enter the interaction strength: "))
    Tmin=float(input("Enter the minimum temperature: "))
    Tmax=float(input("Enter the maximum temperature: "))
    dT=float(input("Enter the temperature step: "))
    print(100*"=")
    return magnetisationTemperature(N,J,Tmin,Tmax,dT)

#magnetisationTemperature(10,1,1.5,3.2,0.01)
terminalMagnetisationTemp()