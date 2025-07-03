# Ising Model Solver
This project provides the tools to various Ising models based on the Hamiltonian

$$\hat{H}=-J\sum_{\langle i,j\rangle} \hat{S}_z^i\hat{S}_z^j-h\sum_i\hat{S}_z^i$$

Included Models and Methods:
- The 1D Quantum Ising model, solved using exact diagonalisation
- The 2D Ising model, solved using a classical Metropolis Monte-Carlo approximation

The pdf in the "Ising Model Notes" folder contain more detailed theoretical background information on the subject matter. Preset solutions are available for viewing in the "Saved Runs" folder. 

## 1D Ising Model 
Run 1DIsing.py to solve for the 1D Ising model under some external magnetic field. 

## 2D Ising Model 
Run twoDIsing_MagField.py to solve the 2D Ising model under some external magnetic field. Run tempMagnetisation.py to plot temperature against magnetisation for the lattice. 

## Notes and Saved Runs
To view the fully rendered notes, visit the "Ising Model Notes" Folder and open "IsingNotes.pdf". Pre-saved solutions/quenches to the 2D Ising model are available as videos in the "Saved Runs Folder" and are too large
to be previewed on github. 
