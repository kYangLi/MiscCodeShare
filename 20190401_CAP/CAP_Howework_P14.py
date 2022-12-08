#!/usr/bin/python3
#
# Date: 2018.12.25
# Author: Liyang@Tsinghua
# Descripution:
# This python3 code is designed to solve the homework of Cold Atom Physics 
#   problem 14

import sys 
import numpy as np
import matplotlib.pyplot as plt

###############################
### Sub Function Defination ###
###############################
def solute_double_well_H(partical_quantity, U, t):
  ## Initialize the Hamiltonian 
  H = np.zeros(shape=(partical_quantity+1,partical_quantity+1),dtype=float)
  ## Go through every state for this system
  # Tips: Since the Hamiltonlian perserve the quantity of paricles, we could 
  #       use the number of particles in one well to represent the whole state! 
  #
  #       n1:   0 1 2 3 ... N
  #           -               -
  #           | * * * * ... * |
  #           | * * * * ... * |
  #       H = | .           . |
  #           | .           . |
  #           | * * * * ... * |
  #           -               -
  for particle_well_1 in range(partical_quantity+1):
    # Off-diagonal term 
    # H_{n1,n1-1}
    if (particle_well_1 - 1 >= 0):
      H[particle_well_1,particle_well_1-1] = -1.0 * t * \
        (particle_well_1*(partical_quantity-particle_well_1+1))**(0.5)
    # H_{n1,n1+1}
    if (particle_well_1 + 1 <= partical_quantity):
      H[particle_well_1,particle_well_1+1] = -1.0 * t * \
        ((particle_well_1+1)*(partical_quantity-particle_well_1))**(0.5)
    # Diagonal term 
    # Tips: Let n1 + n2 = N
    #       n1(n1-1) + n2(n2-1)
    #       = n1(n1-1) + (N-n1)(N-n1-1)
    #       = N(N-1) - 2n1(N-n1)
    # Define some various to simplify and short the expression
    N = partical_quantity
    n1 = particle_well_1
    H[particle_well_1,particle_well_1] = (U/2.0) * (N*(N-1) - 2*n1*(N-n1))  
  eigenvalue_set, eigenvector_set = np.linalg.eig(H)
  eigenvector_set = np.transpose(eigenvector_set)
  ## Sort the eigenvalue and eigenvector
  eigenvalue_sorted_index_list = np.argsort(eigenvalue_set)
  sorted_eigenvalue_set = np.zeros(partical_quantity+1)
  sorted_eigenvector_set = np.zeros((partical_quantity+1,partical_quantity+1))
  for eigenvalue_index in range(partical_quantity+1): 
    sorted_eigenvalue_set[eigenvalue_index] = \
      eigenvalue_set[eigenvalue_sorted_index_list[eigenvalue_index]]
    sorted_eigenvector_set[eigenvalue_index] = \
      eigenvector_set[eigenvalue_sorted_index_list[eigenvalue_index]]
  ## Return the result
  return sorted_eigenvalue_set, sorted_eigenvector_set

def cal_wavefunction_distribution(wave_vector):
  # About the defination of l, please check the textbook P100, Eq(7.9).
  basic_vector_quantity = np.size(wave_vector)
  list_of_l = np.zeros(basic_vector_quantity,dtype=int)
  distribution_of_wavefunction = np.zeros(basic_vector_quantity)
  for n1_index in range(basic_vector_quantity):
    list_of_l[n1_index] = n1_index - ((basic_vector_quantity-1)/2)
    distribution_of_wavefunction[n1_index] = (wave_vector[n1_index])**2

  return list_of_l, distribution_of_wavefunction

def cal_density_matrix(wave_vector):
  basic_vector_quantity = np.size(wave_vector)
  density_11 = 0
  density_12 = 0
  density_22 = 0
  for n1_index in range(basic_vector_quantity):
    density_11 += ((wave_vector[n1_index])**2) * n1_index 
    density_22 += ((wave_vector[n1_index])**2) * \
                   (basic_vector_quantity-1-n1_index) 
    if n1_index < basic_vector_quantity-1:
      density_12 +=  wave_vector[n1_index] * wave_vector[n1_index+1] * \
                     ((n1_index+1)*(basic_vector_quantity-1-n1_index))**0.5
  density_21 = density_12
  density_matrix = np.array([[density_11, density_12],
                             [density_21, density_22]])
  return density_matrix

def cal_relative_fluctuations(wave_vector):
  basic_vector_quantity = np.size(wave_vector)
  average_of_squre = 0
  average = 0
  for n1_index in range(basic_vector_quantity):
    average_of_squre += ((wave_vector[n1_index]**2) * \
                         ((basic_vector_quantity-1-2*n1_index)**2))
    average += ((wave_vector[n1_index]**2) * \
                (basic_vector_quantity-1-2*n1_index))
  relative_fluctuations = average_of_squre - average**2

  return relative_fluctuations

####################
### Main Process ###
####################
if __name__ == "__main__":  
  ### Perpare Before Start ###
  ## Python version check 
  python_version = sys.version
  if not ("3.6" in python_version):
    print("Python3.6 Supported Only!")
    exit("Exit...") 

  ### Parameter read in ###
  print(">>> Parameters Read In <<<")
  print("Please input the particle number:")
  particle_quantity = (int)(input('> '))
  print("Please input the value of interaction energy U:")
  U = (float)(input('> '))
  print("Please input the value of overlap t:")
  t = (float)(input('> '))
  print("Please input the min particle number for gap calculation:")
  min_particle_quantity = (int)(input('> '))
  print("Please input the max particle number for gap calculation:")
  max_particle_quantity = (int)(input('> '))
  print("Please input the particle number change step for gap calculation:")
  particle_quantity_step = (int)(input('> '))
  
  ### Main calculation ###
  print(" ")
  print(">>> Main Calculation <<<")
  print(">> Calculating Question 1...")
  ## Question 1: WF, DM, RPF
  # Get the Eigenstate and Eigenvalue 
  eigenvalue_set, eigenvector_set = \
    solute_double_well_H(particle_quantity, U, t)
  # Get the ground state
  ground_state_vector = np.array(eigenvector_set[0])
  ground_state_energy = eigenvalue_set[0]
  for eigen_index in range(1, np.size(eigenvalue_set)+1):
    # If the difference in energy between two state is less than 0.001(U-t),
    #   then we regard it as the degeneracy state.
    if np.abs((eigenvalue_set[eigen_index]-eigenvalue_set[0])/(U-t)) > 0.001:
      break 
    ground_state_vector += eigenvector_set[eigen_index]
  # WF: Wave Function
  list_l, distribution_of_wavefunction = \
    cal_wavefunction_distribution(ground_state_vector) 
  # DM: Density Matrix
  density_matrix = cal_density_matrix(ground_state_vector)
  temp_eigenvalue, temp_eigenvector = np.linalg.eig(density_matrix)
  diag_density_matrix = \
    np.array([[temp_eigenvalue[0],0],[0,temp_eigenvalue[1]]]) 
  # RPF: Relative Particle Fluctuations
  relative_fluctuations = cal_relative_fluctuations(eigenvector_set[0])

  print(">> Calculating Question 2...") 
  ## Question 2: Gap
  particle_quantity_list = np.array([]) 
  energy_gap_list = np.array([]) 
  for particle_quantity_for_gap in range(min_particle_quantity,
                                         max_particle_quantity+1,
                                         particle_quantity_step):
    # Precess print
    process_print = '> Particle Quantity: ' + str(particle_quantity_for_gap)
    delete_print = '\b' * len(process_print)
    sys.stdout.write(process_print)
    sys.stdout.flush()
    sys.stdout.write(delete_print)
    # Gap calculate
    particle_quantity_list = np.append(particle_quantity_list,
                                       particle_quantity_for_gap)
    eigenvalue_set, eigenvector_set = \
      solute_double_well_H(particle_quantity_for_gap, U, t)
    ground_state_energy = eigenvalue_set[0]
    for eigen_index in range(1, np.size(eigenvalue_set)):
      # If the difference in energy between two state is less than 0.001,
      #   then we regard it as the degeneracy state.
      excited_state_energy_index = eigen_index
      if np.abs((eigenvalue_set[eigen_index]-eigenvalue_set[0])/(U-t)) > 0.001:
        break
    excited_state_energy = eigenvalue_set[excited_state_energy_index]
    energy_gap = excited_state_energy - ground_state_energy
    energy_gap_list = np.append(energy_gap_list, energy_gap)

  ### Result Output ###
  print("")
  print("")
  print(">>> Output Result <<<")
  print('-------------------------------------')
  print('=> Density Matrix:')
  print(diag_density_matrix)
  print("")
  print("=> Relative Fluctuations:")
  print(relative_fluctuations)
  print("")
  #
  file_name = "Psi-l-U" + str(U) + "t" + str(t) + \
              "N" + str(particle_quantity) + ".jpg"
  print("=> Plotting the ", file_name)
  plt.plot(list_l, distribution_of_wavefunction)
  plt.xlabel("l")
  plt.ylabel("$|\Psi_l|^2$")
  plt.grid()
  plt.savefig(file_name)
  plt.show()
  #
  file_name = "GE-N-U" + str(U) + "t" + str(t) + ".jpg"
  print("=> Plotting the ", file_name)
  plt.plot(particle_quantity_list, energy_gap_list)
  plt.xlabel("Particle Quantity")
  plt.ylabel("Energy Gap")
  plt.grid()
  plt.savefig(file_name)
  plt.show()
  print('-------------------------------------')
