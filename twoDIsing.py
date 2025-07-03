import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import imageio
import random as rd
import os

def gridGen(N:int)->list[list[int]]:
    """
    Generate a 2D lattice of random spins. Use dimensionless spin +/-1 for up/down.  
    """
    i=0 
    j=0 
    outputGrid=[]
    while i<N:
        outputGrid.append([])
        j=0
        while j<N: 
            outputGrid[i].append(rd.choice([-1,1]))
            j+=1
        i+=1
    return outputGrid

def hamiltonian(N:int,J:float,lattice:list[list[int]] | None =None):
    """
    Return the Ising energy of an N x N lattice with coupling strength J. 
    Apply periodic boundary conditions (loop back to i=0 and j=0 when at edge)
    """
    if lattice==None:
        lattice=gridGen(N)
    i=0 
    j=0 
    totalEnergy=0
    while i<N:
        j=0
        while j<N:
            totalEnergy+=lattice[i][j]*lattice[i][(j+1) % N]
            totalEnergy+=lattice[i][j]*lattice[(i+1) % N][j]
        #    print("(", i, ",",j,")")
            j+=1
        i+=1
    return -J*totalEnergy
def magnetisation(N:int,lattice:list[list[int]] | None =None):
    """
    Return the magnetisation of an NxN Ising lattice 
    """
    if lattice is None:
        lattice=gridGen(N)
    mag=0 
    i=0 
    j=0 
    while i<N:
        j=0 
        while j<N:
            mag+=lattice[i][j]
            j+=1
        i+=1
    return mag

def chars(N,J):
    lattice=gridGen(N)
    energy=hamiltonian(N,J,lattice)
    mag=magnetisation(N,lattice)
    return energy, mag

def spinPlotter(lattice):
    i=0 
    j=0 
    N=len(lattice)
    while i<N:
        j=0
        while j<N:
            if lattice[i][j]<0:
                plt.plot(i,j,".",color="red")
            else:
                plt.plot(i,j,".",color="blue")
            j+=1
        i+=1
    plt.show()

def lattice_to_image(lattice, size=4):
    fig = Figure(figsize=(size, size))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.imshow(lattice, cmap='bwr', vmin=-1, vmax=1)
    ax.axis('off')
    canvas.draw()
    img = np.frombuffer(canvas.tostring_rgb(), dtype='uint8')
    img = img.reshape(canvas.get_width_height()[::-1] + (3,))
    plt.close(fig)
    return img

def metropolisFlips(N:int,J:float, T:float, maxIters:int=10000, fpsCustom:int=100, verbosity:bool=True, lattice:list[list[int]] | None =None):
    """
    Return the lattice in its thermal equillibrium state under temperature T and coupling strength J. 
    Solve using the Monte Carlo method, implementing the Metropolis model of interaction* 
    
    *: A spin flip is guaranteed to happen when it causes the overall energy to be lower. It can also happen if it raises the energy with a probabilty of P=exp(dE/T) where dE is the 
    energy raised by the flip

    Note: kB=1, express T in units of the coupling strength J
    """
    frames=[]
    energies=[]
    iterationList=[]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    videoPath = os.path.join(script_dir, "isingSim_2D.mp4")
    if lattice==None:
        lattice=gridGen(N)
    energyStable=0 
    iters=0
    while energyStable<=1000 and iters<=maxIters:
        i=rd.randint(0,N-1)
        j=rd.randint(0,N-1)
        if verbosity:
            print("Trying to flip ", (j,i))
        newEnergy= hamiltonian(N,J,lattice)
        energies.append(newEnergy)
        iterationList.append(iters)
        f = (lattice[(i+1)%N][j] +lattice[(i-1)%N][j] +lattice[i][(j+1)%N] +lattice[i][(j-1)%N])
        dE=2*J*lattice[i][j]*f
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
            energyStable=0
            lattice[i][j]=-lattice[i][j]
            if verbosity:
                print("Flip successful")
        if verbosity:
            for line in lattice:
                print(line)
        print(30*"\n","Solver running. Do not close or interrupt program.\n Iteration ", iters, "out of a maximum possibility of ", maxIters,"\n Equilibrium not reached yet. Current energy=", newEnergy)
        frames.append(lattice_to_image(lattice))
        if verbosity:
            print("Energy change from flip: ", dE)
            print("Energy: ", newEnergy)
            print("No Change Counter: ", energyStable)
            print("Frame added, current frame count: ",len(frames) )
   #     print(30*"=")
        iters+=1
    print(30*"\n")
    print(100*"-")
    if iters>=maxIters:
        print("Solving Complete. Equilibrium state not reached before maximum iteration limit. ")
    else: 
        print("Solving Complete.")
        if T<0.0001:
            print("Equilibrium state reached. Found ground state energy eigenvalue E=", newEnergy)
        else:
            print("Equilibrium state reached at an energy of E=", newEnergy)
        print("Iterations before convergence: ", iters)
    print(100*"-")
    print("\nPost-processing solution:")
    print("Constructing video representation...")
    print("Plotting iteration-by-iteration energies...")
    imageio.mimsave(videoPath, frames, fps=fpsCustom)
    print("Video saved to ", videoPath)
   # spinPlotter(lattice)
    plt.plot(iterationList,energies,color="black")
    plt.xlabel("Iterations")
    plt.ylabel("Energy of Lattice")
    print("Post-processing complete. See output plot for energy evolution and output video for lattice video. ")
    plt.show()
    return lattice

def terminalMetropolis():
    print(100*"=")
    print("Welcome to the Ising model solver")
    preview=int(input("Enter 1 to view a preset simulation or any other key to continue: "))
    print(100*"=")
    if preview==1:
        N=50
        J=1
        T=0.0000000000001
        maxIters=30000
        fpsCustom=5000
        fancy=False
    else:
        print("Simulation Conditions: ")
        N=int(input("Enter the size of the lattice (Number of particles on each side): "))
        J=float(input("Enter the coupling constant of the particles: "))
        T=float(input("Enter the temperature (In multiples of the coupling constant): "))
        print(100*"-")
        print("Computational Parameters: ")
        maxIters=int(input("Enter the maximum allowed iterations (Recommended 60000): "))
        fpsCustom=int(input("Enter the fps of the output video (Recommended 100): "))
        fancy=False 
        debug=int(input("Enter 1 to enable debug mode or any other key to continue: "))
        if debug==1:
            fancy=True 
    return metropolisFlips(N,J,T,maxIters,fpsCustom, fancy)

#metropolisFlips(10,1,0.0000000001,30000,5000,False)
#terminalMetropolis()   