
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import imageio
import random as rd
import os
import twoDIsing_NoMagField as tDI

#Git push test

def metropolisMagnetisation(N:int,J:float, T:float, H:float, maxIters:int=50000, fpsCustom=100, verbosity:bool=True, lattice:list[list[int]] | None =None):
    """
    Return the lattice in its thermal equillibrium state under temperature T and coupling strength J. 
    Solve using the Monte Carlo method, implementing the Metropolis model of interaction* 
    
    *: A spin flip is guaranteed to happen when it causes the overall energy to be lower. It can also happen if it raises the energy with a probabilty of P=exp(dE/T) where dE is the 
    energy raised by the flip

    Note: kB=1, express T in units of the coupling strength J

    H is the external magnetic field strength
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    videoPath = os.path.join(script_dir, "isingSim_2D.mp4")
    energies=[]
    frames=[]
    iterationList=[]
    magnetisations=[]
    if lattice==None:
        lattice=tDI.gridGen(N)
    energyStable=0 
    magnetisation=0
    iters=0
    while energyStable<=69*N and iters<=maxIters:
        i=rd.randint(0,N-1)
        j=rd.randint(0,N-1)
        if verbosity:
            print("Trying to flip ", (j,i))
        iterationList.append(iters)
        f = (lattice[(i+1)%N][j] +lattice[(i-1)%N][j] +lattice[i][(j+1)%N] +lattice[i][(j-1)%N])
        dE= 2*lattice[i][j]*(J*f+H) 
        if verbosity:
            print("Energy change if successsful:", dE)
        if dE>0:
            b=1/T
            transProb=np.exp(-b*dE)
            r=rd.random()
            if r>transProb:
                dE=0
                energyStable+=1
                if verbosity:
                    print("Flip failed")
            else:
                energyStable=0
                lattice[i][j]=-lattice[i][j]
                if verbosity:
                    print("Flip successful")
        elif dE<=0:
            lattice[i][j]=-lattice[i][j]
            energyStable=0
            if verbosity:
                print("Flip successful")
        if verbosity:
            for line in lattice:
                print(line)
        newEnergy= tDI.hamiltonian(N,J,lattice)
        energies.append(newEnergy)
        print(30*"\n","Solver running for N=",N,". Do not close or interrupt program.\n Iteration ", iters, "out of a maximum possibility of ", maxIters,"\n Equilibrium not reached yet.\n Current status: Energy=", newEnergy, "|| Magnetisation=", magnetisation)
        frames.append(tDI.lattice_to_image(lattice))
        if verbosity:
            print("Energy change from flip: ", dE)
            print("Energy: ", newEnergy)
            print("No Change Counter: ", energyStable)
   #     print(30*"=")
        magnetisation=tDI.magnetisation(N,lattice)
        magnetisations.append(magnetisation)
        iters+=1
    print(30*"\n")
    print(160*"-")
    if iters>=maxIters:
        print("Solving Complete. Equilibrium state not reached before maximum iteration limit. Displaying conditions of the last iteration.")
        print("Energy: E=", newEnergy, " || magnetisation: M=", magnetisation)
    else: 
        print("Solving Complete. Equilibrium state reached at the conditions shown below.")
        print("Energy: E=", newEnergy, " || Magnetisation: M=", magnetisation)
        print("Iterations before convergence: ", iters)
    print(160*"-")
    print("\nPost-processing solution:")
    print("Constructing video representation...")
    print("Plotting iteration-by-iteration energies...")
    imageio.mimsave(videoPath, frames, fps=fpsCustom)
    print("Video saved to ", videoPath)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))
    ax1.plot(iterationList,magnetisations,color="black")
    ax1.set_xlabel("Iterations")
    ax1.set_ylabel("Magnetisation")
    ax2.plot(iterationList,energies,color="black")
    ax2.set_xlabel("Iterations")
    ax2.set_ylabel("Energy")
    fig.suptitle("Energies and Magnetisation vs Iterations")
    print("Post-processing complete. See output plot for energy evolution and output video for lattice video. ")
    i=0 
    j=0
    while i<N:
        j=0
        while j<N:
            if lattice[i][j]<0:
                plt.plot(j,i,".",color="red")
            else:
                plt.plot(j,i,".",color="blue")
            j+=1
        i+=1
    plt.show()
    return magnetisation

def terminalMetropolis():
    print(100*"=")
    print("Welcome to the Ising model solver")
    preview=int(input("Enter 1 to view a preset simulation or any other key to continue: "))
    print(100*"=")
    if preview==1:
        N=50
        J=1
        H=50
        T=0.0000000000001
        maxIters=30000
        fpsCustom=5000
        fancy=False
    else:
        print("Simulation Conditions: ")
        N=int(input("Enter the size of the lattice (Number of particles on each side): "))
        J=float(input("Enter the coupling constant of the particles: "))
        T=float(input("Enter the temperature (In multiples of the coupling constant): "))
        H=float(input("Enter the strength of the external magnetic field: "))
        print(100*"-")
        print("Computational Parameters: ")
        maxIters=int(input("Enter the maximum allowed iterations (Recommended 86000): "))
        fpsCustom=int(input("Enter the fps of the output video (Recommended 100): "))
        fancy=False 
        debug=int(input("Enter 1 to enable debug mode or any other key to continue: "))
        if debug==1:
            fancy=True 
    return metropolisMagnetisation(N,J,T,H,maxIters,fpsCustom, fancy)


#High magnetic field magnitude => converge (all polarised in same direction) faster, + => up, -=> down, |M| converges to N^2 
#Magnitisation magnetic field near 0 => more chaotic, less neatly polarised
terminalMetropolis()