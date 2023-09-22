import numpy as np
import os
import pandas as pd
import sys
from tkinter import filedialog
import tkinter as tk
from tkinter import messagebox
import warnings  # Import the warnings module

#------------ Input Parameters for EAL-:

rho = 2.5  # density in g/cm^3
N_v = 16  # valence electrons per molecule
MZ = 80  # molecular mass
E_gap = 3.2  # band gap in eV

#------------ File Options ------------:


#file in:
root = tk.Tk()
root.withdraw()
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

with open(f_in, 'r', encoding='latin-1') as file:
    lines = file.readlines()
    lines = [line.strip() for line in lines if line.strip()]


df = pd.DataFrame([line.split() for line in lines])


in_array = df.to_numpy()

# Read the file using pandas
#df = pd.read_csv(f_in, sep=' ', header=0, encoding='latin-1')
#df = pd.read_csv(f_in, sep=',', header=0, encoding='latin-1')
df = pd.read_csv(f_in, sep='\t', header=0, encoding='latin-1')




# Convert DataFrame to NumPy array
in_array = df.to_numpy()

#get # of energies and species:
n_hv = np.size(in_array, axis=0) - 1
n_spec = np.size(in_array, axis=1) - 1

#get hv and species (BE):
hv = []
for k in range(0, n_hv):
    hv.append(in_array[k + 1, 0])

spec = []
for j in range(0, n_spec):
    spec.append(in_array[0, j + 1])

#create measured array - kxj array:
Meas = np.empty((n_hv, n_spec))
for k in range(0, n_hv):
    for j in range(0, n_spec):
        Meas[k, j] = in_array[k + 1, j + 1]
        
#------------ Parameters ------------:

# Define a function to get user inputs through a pop-up window
def get_user_inputs():
    root = tk.Tk()
    root.withdraw()

    # Ask for file input
    f_in = filedialog.askopenfilename()

    # Ask for other parameters
    theta = float(input("Enter theta: "))
    delta = float(input("Enter delta value: "))
    max_it = int(input("Enter maximum number of iterations: "))
    max_a_it = int(input("Enter maximum number of alpha iterations: "))
    conv_crit = float(input("Enter convergence criterion: "))
    n_smooth = int(input("Enter number of smoothing iterations: "))
    stepsize = float(input("Enter step size (in nm): "))

    return f_in, theta, delta, max_it, max_a_it, conv_crit, n_smooth, stepsize


#------------ Functions ------------:
f_in, theta, delta, max_it, max_a_it, conv_crit, n_smooth, stepsize = get_user_inputs()

def calcEAL(theta, Meas, MZ, N_model, E_gap):
    # calculate EAL using reference 4
    N_hv, N_spec = np.shape(Meas)
    EAL = np.zeros((N_hv, N_spec))  # Initialize EAL array

    for k in range(N_hv):
        for j in range(N_spec):
            term1 = MZ * calc[k, j]  # Replace calc with Meas
            term2 = N_model[k, j] * E_gap / (N_v * rho)  # Corrected N_model indexing
            term3 = (1 - np.exp(-theta * Meas[k, j])) / (1 - np.exp(-theta * calc[k, j]))  # Replace calc with Meas
            EAL[k, j] = term1 * term2 * term3

    return EAL



"""
def calcEAL(theta, Meas, MZ, N_model, E_gap):    
   # calculate EAL using reference 4
   # N_hv, N_spec = np.shape(Meas)
   # EAL = np.zeros((N_hv, N_spec))
    
    #for k in range(N_hv):
    #    for j in range(N_spec):
    #        term1 = MZ * calc[k, j]
    #       term2 = N_model * E_gap / (N_v * rho)
    #        term3 = (1 - np.exp(-theta * Meas[k, j])) / (1 - np.exp(-theta * calc[k, j]))
    #        EAL[k, j] = term1 * term2 * term3
    
    #calc = np.zeros_like(Meas)
    calc = np.ones((N_hv, N_spec))  # Initialize calc with the same shape as Meas
    M_j = np.zeros_like(Meas)
    #N_model = np.zeros_like(Meas)
    N_model = np.ones((N_hv, N_spec))  # Corrected initialization of N_model
    alpha = 1.0
    chisq_old = 1e30
    entropy_old = 1e30

    # Iterative algorithm
    for it in range(max_it):
        for a_it in range(max_a_it):
            #calc_new = calcModel(theta, MZ, N_model, rho, E_gap)
            calc_new = calcCalc(theta, calc, M_j, N_model)  # Corrected function name
            chisq = calcChisq(Meas, calc_new, theta)

            M_j_new = M_j - alpha * (calc_new - calc) * calc_new / (calc_new ** 2 + delta)
            N_model_new = N_model - alpha * (calc_new - calc) * M_j_new / (calc_new ** 2 + delta)

            entropy_new = calcEntropy(Meas, calc_new)

            if entropy_new > entropy_old:
                alpha /= 10
                if alpha < 1e-6:
                    alpha = 1e-6
            else:
                alpha *= 10
                if alpha > 1e6:
                    alpha = 1e6

            if entropy_new < entropy_old:
                calc = calc_new
                M_j = M_j_new
                N_model = N_model_new
                entropy_old = entropy_new
                chisq_old = chisq
    
    return EAL
"""
###new func
def calcModel(theta, MZ, N_model, rho, E_gap):
    # calculate calc using reference 1
    N_hv, N_spec = np.shape(N_model)
    calc_new = np.zeros((N_hv, N_spec))

    for k in range(N_hv):
        for j in range(N_spec):
            term1 = 1 - np.exp(-theta * MZ)
            term2 = (1 - term1) / (1 - np.exp(-theta * N_model[0, j]))
            calc_new[k, j] = np.sum(rho * N_v * term2)

    return calc_new
    
def calcCalc(theta, calc, M_j, N_model):
    # calculate calc using reference 1
    N_hv, N_spec = np.shape(calc)
    calc_new = np.zeros((N_hv, N_spec))
    
    for k in range(N_hv):
        for j in range(N_spec):
            term1 = 1 - np.exp(-theta * calc[k, j])
            term2 = (1 - term1) / (1 - np.exp(-theta * M_j[k, :]))
            calc_new[k, j] = np.sum(N_model * term2)
    
    return calc_new

def calcChisq(Meas, calc, theta):
    # calculate chi-squared value
    N_hv, N_spec = np.shape(Meas)
    chisq = 0
    
    for k in range(N_hv):  # Now iterate over the correct axis
        for j in range(N_spec):  # Now iterate over the correct axis
            term1 = Meas[k, j] - calc[k, j]
            if np.isclose(calc[k, j], 0):
                term2 = 0
            else:
                term2 = term1 ** 2 / (1 - np.exp(-theta * calc[k, j]) + 1e-8)
            chisq += term2            
    
    return chisq

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



#------------ Main Code ------------:

# initialize parameters
N_hv, N_spec = np.shape(Meas)
chisq_old = 1e20
calc = np.ones((N_hv, N_spec))
M_j = np.ones((N_hv, 1))
N_model = np.ones((1, N_spec))
alpha = 1
entropy_old = calcEntropy(Meas, calc)

# Iterative algorithm
for it in range(max_it):
    for a_it in range(max_a_it):
        #uncommented for line below -- calc_new = calcModel(theta, MZ, N_model, rho, E_gap)
        calc_new = calcCalc(theta, calc, M_j, N_model)
        chisq = calcChisq(Meas, calc_new, theta)
        
        if chisq > chisq_old:
            alpha /= 10
            if alpha < 1e-6:
                alpha = 1e-6
        elif chisq < chisq_old:
            alpha *= 10
            if alpha > 1e6:
                alpha = 1e6
        
        M_j_new = M_j - alpha * (calc_new - calc) * calc_new / (calc_new ** 2 + delta)
        N_model_new = N_model - alpha * (calc_new - calc) * M_j_new / (calc_new ** 2 + delta)

        entropy_new = calcEntropy(Meas, calc_new)

        if entropy_new > entropy_old:
            alpha /= 10
            if alpha < 1e-6:
                alpha = 1e-6
        else:
            alpha *= 10
            if alpha > 1e6:
                alpha = 1e6

        if entropy_new < entropy_old:
            calc = calc_new
            M_j = M_j_new
            N_model = N_model_new
            entropy_old = entropy_new
            chisq_old = chisq
    
    EAL = calcEAL(theta, Meas, MZ, N_model, E_gap)
    chisq_final = calcChisq(Meas, calc_new, theta)
    entropy_final = calcEntropy(Meas, calc_new)   
    
    def display_results(EAL, chisq_final, entropy_final, hv, spec):
        msg = f"Final Chisq: {chisq_final}\nFinal Entropy: {entropy_final}"
        for k, hv_value in enumerate(hv):
            msg += f"\n\nEAL at hv={hv_value}:\n"
            for j, spec_value in enumerate(spec):
                msg += f"{spec_value}: {EAL[k][j]}\n"
                
            messagebox.showinfo("Results", msg)

    # Call the function after the calculations are done
    display_results(EAL, chisq_final, entropy_final, hv, spec)

    
    """""
    # output results
    #test outputs to test code 
    print("theta, Meas, MZ, N_model, E_gap", theta, Meas, MZ, E_gap)
    
    print("Final Results:")
    print("Final Chisq: ", chisq_final)
    print("FInal Entropy: ", entropy_final)
    """

#return not in function? 
# pilace in calc eal then call function towards end of prog?
#return EAL, chisq_final, entropy_final

