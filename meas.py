import numpy as np
import os
import pandas as pd
import sys
from tkinter import filedialog
import tkinter as tk
from tkinter import messagebox
import warnings  # Import the warnings module

#------------ Input Parameters for EAL-:


#make input in future
rho = 2.5  # density in g/cm^3
N_v = 16  # valence electrons per molecule
MZ = 80  # molecular mass
E_gap = 3.2  # band gap in eV

#------------ File Options ------------:

#file in:
root = tk.Tk()
root.withdraw()
messagebox.showinfo("Prompt", "Open Test-Data")
f_in = filedialog.askopenfilename()

#file prefix:
f_pref = os.path.splitext(f_in)[0]

#out files:
f_par = f_pref + '.par'   # fitting parameters
f_alpha = f_pref + '.alf'
f_calc = f_pref + '_calc.vke'    # calculated VKE prof.
f_prof = f_pref + '.dep'   # final depth profile
f_eal = f_pref + '.eal'  # EALs

#------------ Input Arrays ------------:

    #this input sections should recongize theta and apply it to prog


with open(f_in, 'r', encoding='latin-1') as file:
    lines = file.readlines()
    lines = [line.strip() for line in lines if line.strip()]


df = pd.DataFrame([line.split() for line in lines])


in_array = df.to_numpy()

# Read the file using pandas
df = pd.read_csv(f_in, sep='\t', header=0, encoding='latin-1')


# Convert DataFrame to NumPy array
in_array = df.to_numpy()

#get # of energies and species:
N_hv = np.size(in_array, axis=0) - 1
n_spec = np.size(in_array, axis=1) - 1

#------------ End Of Input Arrays ------------:

#get hv and species (BE):
hv = []
for k in range(0, N_hv):
    hv.append(in_array[k + 1, 0])

spec = []
for j in range(0, n_spec):
    spec.append(in_array[0, j + 1])

#create measured array - kxj array:
Meas = np.empty((N_hv, n_spec))
for k in range(0, N_hv): #last point for capat
    for j in range(0, n_spec):
        Meas[k, j] = in_array[k + 1, j + 1]

# Print out Meas array
print("Print Meas array:")
for k in range(N_hv):
    for j in range(n_spec):
        print(f"{Meas[k, j]:.4f}", end="\t")
    print()
   
   
     
#print M_J
print("Printing MJ")
M_j = np.ones((N_hv, 1))
print(M_j)


def calcEntropy(Meas, calc):
    # calculate entropy value
    N_hv, N_spec = np.shape(Meas)
    entropy = 0
    
    for k in range(N_hv):
        for j in range(N_spec):
            term1 = Meas[k, j] * np.log(Meas[k, j])
            if np.isclose(calc[k, j], 0):
                term2 = 0
            else:
                term2 = calc[k, j] * np.log(calc[k, j])
            term3 = Meas[k, j] - calc[k, j]
            entropy += term1 - term2 - term3
    
    return entropy


calc_new = np.zeros((N_hv, n_spec))


N_hv, N_spec = np.shape(Meas)
chisq_old = 1e20
calc = np.ones((N_hv, N_spec))
M_j = np.ones((N_hv, 1))
N_model = np.ones((1, N_spec))
alpha = 1
entropy_old = calcEntropy(Meas, calc)

delta = 15

print("pritning mj new")
M_j_new = M_j - alpha * (calc_new - calc) * calc_new / (calc_new ** 2 + delta)

print(M_j_new)

