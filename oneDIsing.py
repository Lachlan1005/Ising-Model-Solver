import numpy as np
import matplotlib.pyplot as plt

Sx = 0.5 * np.array([[0, 1], [1, 0]])  # Pauli X
Sy = 0.5 * np.array([[0, -1j], [1j, 0]])  # Pauli Y
Sz = 0.5 * np.array([[1, 0], [0, -1]])  # Pauli Z
I=np.eye(2)

def qubits(N:int,dir:int):
    """
    Initialise the spins of N qubits arranged in a 1D straight line lattice. 
    Return a list with N items, representing the spins of each qubit. 
    """
    spins=[]
    n=0 
    while n<N:
        prod=[]
        i=0 
        while i<N:
            if i==n:
                if dir==1:
                    prod.append(Sx)
                elif dir==2:
                    prod.append(Sy)
                elif dir==3:
                    prod.append(Sz)
            else:
                prod.append(I)
            i+=1
        j=1
        spin=prod[0]
        while j<N:
            spin=np.kron(spin,prod[j])
            j+=1
        spins.append(spin)
        n+=1
    return spins 

def qubitInteractions(N:int,J:float):
    """
    Return the qubit interaction term of the Ising hamiltonian for N qubits and coupling constant J
    """
    H0=0
    n=0 
    spins=qubits(N,3)
    while n<N-1:
        H0+=spins[n] @ spins[n+1]
        n+=1
    return -J*H0 


def fieldInteractions(N:int,h:float):
    """
    Return the external magnetic field interaction term of the Ising hamiltonian for N qubits and field strength h
    """
    H0=0 
    n=0
    spins=qubits(N,1)
    while n<N:
        H0+=spins[n]
        n+=1
    return -h*H0 

def hamiltonian(N:int,J:float,h:float)->dict[any]:
    """
    Return a dict with the associated {"Hamiltonian": "" , Energy Levels: "" } for a 1D ising lattice with N qubits with coupling constant J under 
    a magnetic field with strength h
    """
    H=qubitInteractions(N,J)+fieldInteractions(N,h)
    E=np.linalg.eigvalsh(H)
    sol={"Your input [N,J,h]":[N,J,h], "Hamiltonian":H, "Energy Eigenvalues":E}
    for key in sol:
        print("\n\n",key,"\n" ,sol[key], "\n=======================")
    return sol

def groundState(N:int, J:float, h:float, Jorh:int, max:float, step:float):
    """
    Plots the gound state energy levels against either J or h from the input variable to max, incrementing with step size specified by the input parameter step
    if Jorh==0 -> plot for J
    if Jorh==1 -> plot for h
    """
    x=[]
    E=[]
    J0=J
    h0=h
    i=1
    if Jorh==0:
        totali=(max-J)/step
        var="Qubit Interaction Coupling Parameter"
        while J<=max:
            H=hamiltonian(N,J,h)
            x.append(J)
            E.append(H["Energy Eigenvalues"][0])
            #plt.plot(x,E,"-",color="black")
            print("Progress:", (i*100)/(totali), "%")
            i+=1
            J+=step
    if Jorh==1:
        var="External Magnetic Field Strength"
        totali=(max-h)/step
        while h<=max:
            H=hamiltonian(N,J,h)
            x.append(h)
            E.append(H["Energy Eigenvalues"][0])
            #plt.plot(x,E,"-",color="black")
            print("Progress:", (i*100)/totali, "%")
            i+=1
            h+=step
    plt.plot(x,E,color="black")
    plt.xlabel(var)
    plt.ylabel("Ground State Energy")
    plt.show()

def spinComp(N:int, J:float, h:float, dir:int):
    H=qubitInteractions(N,J)+fieldInteractions(N,h)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    gStateVal=eigenvalues[0]
    psi=eigenvectors[:, 0]
    n=0
    spinOps=qubits(N,dir)
    expectations=[]
    while n<N: 
        expectation=np.dot(psi.conj().T, spinOps[n] @ psi)
        expectations.append(round(expectation*100000)/100000)
        n+=1
    return expectations

def spin(N:int, J:float, h:float)->list[tuple[float,float]]:
    print("Solver running...")
    qubitList=[]
    n=0 
    spinsX=spinComp(N,J,h,1)
    spinsZ=spinComp(N,J,h,3)
    while n<N:
        qubitList.append((spinsX[n],spinsZ[n]))
        n+=1
    return qubitList

def spinBloch(N,J,h):
    qubitSpins=spin(N,J,h)
    
        
#Vary J
#groundState(11, -1,0.5, 0, 1, 0.01)

#Vary h
#groundState(8, 1,-1, 1, 1, 0.01)

    
print("spin (x,z): ", spin(10,100,1))

