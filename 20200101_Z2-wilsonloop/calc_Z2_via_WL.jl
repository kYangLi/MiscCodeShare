# Author: Yang Li
# Date: 2022.08.31
# Email: lyang.1915@foxmail.com
# Description: Use the HopTB to calculate the wilson loop, and further to calculate the Z2 number. See method in `http://es12.wfu.edu/pdfarchive/oral/ES12_A_Soluyanov.pdf`.

using DelimitedFiles # Support file read and write
using HopTB


### Generate a continuious kpoints path, defined by the begin and end point.
function gen_kpath(k_beg, k_end, kp_num)
  kpath_beg2end = hcat(k_beg, k_end) # Each `column` represents one kpoint.
  return HopTB.Utilities.constructlinekpts(kpath_beg2end, kp_num)
end


### Calculate the zak phase value for each ky point
function calc_wilson_spectrum_for_every_ky!(
  k_path_y, ky_num, kx_end, kx_num, occupy_num, 
  wilson_loop_for_bands, tb_model)
  ##
  for i_ky in 1:ky_num # Index for each ky in ky path.
    println(" +--[sub-do] Point ky #$(i_ky)/$(ky_num)")
    # The KPath in x direction, to integrate calculate one zak phase
    wilson_k_beg = k_path_y[:, i_ky] #[kx(-0.5), ky, kz(0)]
    wilson_k_end = [kx_end, wilson_k_beg[2], 0]#[kx(0.5), ky, kz(0)]
    wilson_kpath_beg2end = hcat(wilson_k_beg, wilson_k_end)
    # All band indies
    all_bands_indices = collect(1:occupy_num)
    # Calculate the zak phase for each ky, all bands once.
    # The `HopTB.Topology.get_wilson_spectrum` fucntion has already sort its output eigvalues of wilson spectra, so no need to sort it here.
    wilson_loop_for_bands[i_ky,:] = 
      HopTB.Topology.get_wilson_spectrum(
        tb_model, all_bands_indices, wilson_kpath_beg2end, kx_num)
  end
  return 0
end


function get_tb_model(data_file, data_type)
  if data_type == "wannier90"
    return HopTB.Interface.createmodelwannier(data_file)
  elseif data_type == "openmx"
    return HopTB.Interface.createmodelopenmx(data_file)
  else
    println("[error] TB model data type must in [wannier90, openmx ]")
    exit()
  end
end

### Calculate the Z2 number from the wilson loop
### - The wilson loop integrate the kx direction.
### - And connect the ky direction to form a loop.
### - The band index are also considered in calculate.
### - Each column is one wilson loop for one band.
function calc_wilson_loop(
  kx_beg, kx_end, kx_num,
  ky_beg, ky_end, ky_num, 
  occupy_num, data_file, data_type)
  ##
  println("[do] Reading TB model...")
  tb_model = get_tb_model(data_file, data_type)
  ##
  println("[do] Calculating the wilson loop for each band and ky-point...")
  # KPath alone the y axis. Later, each kpoint in this path will equip with one zak phase, froming the wilson spectrum. 
  k_path_y = gen_kpath([kx_beg, ky_beg, 0], [kx_beg, ky_end, 0], ky_num)
  # wilson loop for each band, where each column is one wilson loop for one band.
  wilson_loop_for_bands = zeros(ky_num, occupy_num)
  calc_wilson_spectrum_for_every_ky!(
    k_path_y, ky_num, kx_end, kx_num, occupy_num, 
    wilson_loop_for_bands, tb_model)
  ## 
  return k_path_y, wilson_loop_for_bands
end


### First transfer an angle to a complex number, and then vectroize it in the 2D complex space.
### - The angle is in 1 unit (2pi = 3.14...), not in degree.
### - Be care! The final output vector is actually in 3D space not 2D. But latter we will defined a function called `corss_times_z`, it will work even under the two cross vectors not in 3D. So, here the out-of-plane component is not needed.
function vectorize_an_angle(angle)
  exp_iangle = exp(1im * angle)
  return [real(exp_iangle), imag(exp_iangle)]
end


### Calculate the z component of two vectors cross production.
function cross_times_z(vec1, vec2)
  return vec1[1]*vec2[2] - vec1[2]*vec2[1]
end


### Judage the three points' arangment: clockwise (+1) or anti-clockwise(-1).
### - The denotes can be refered in the Slides provided by Xu Runzhang. See also in http://es12.wfu.edu/pdfarchive/oral/ES12_A_Soluyanov.pdf, Page 40.
function get_arangment_chirality(z_m, z_mp1, xbar_n_mp1)
  vec1 = z_m - xbar_n_mp1
  vec2 = z_mp1 - xbar_n_mp1
  return Int(sign(cross_times_z(vec1, vec2)))
end


### Get the middle point of the phases
function get_middle_point(phase1, phase2)
  phase_sml = min(phase1, phase2)
  phase_big = max(phase1, phase2)
  if phase_big - phase_sml >= pi
    middle_phase = (phase_big + phase_sml) / 2
  else
    # The mod(2pi) here is actually not need, since the exp(i*theta) will automatically do that. But still we add it here, to avoid some unnecessary misleading to code readers. In real calculation, pls remove it.
    middle_phase = (((phase_big + phase_sml) / 2) + pi) % 2pi
  end
end


###
function get_necessary_vecs(wilson_loop_for_bands, i_ky, i_band, occupy_num)
  phase_z_m = wilson_loop_for_bands[i_ky, i_band]
  if i_band != occupy_num
    phase_z_mp1 = wilson_loop_for_bands[i_ky, i_band+1]
  else
    phase_z_mp1 = wilson_loop_for_bands[i_ky, 1]
  end
  phase_xbar_n_mp1 = get_middle_point(phase_z_m, phase_z_mp1)
  z_m = vectorize_an_angle(phase_z_m)
  z_mp1 = vectorize_an_angle(phase_z_mp1)
  xbar_n_mp1 = vectorize_an_angle(phase_xbar_n_mp1)
  return z_m, z_mp1, xbar_n_mp1
end


### Sum all ky points to get the value of Z2
function sum_to_Z2(ky_num, occupy_num, wilson_loop_for_bands)
  Z2_number = 0
  bad_point_num = 0
  for i_ky in 1:ky_num # Loop for each ky point
    ky_sign = 1
    for i_band in 1:occupy_num # Loop for each band
      z_m, z_mp1, xbar_n_mp1 = 
        get_necessary_vecs(wilson_loop_for_bands, i_ky, i_band, occupy_num)
      curr_ky_sign = get_arangment_chirality(z_m, z_mp1, xbar_n_mp1)
      ky_sign *= curr_ky_sign
    end
    # bad point skip!
    if ky_sign == 0
      bad_point_num += 1
      continue
    end
    # sum to z2
    ky_sign = div(1-ky_sign, 2) # mapping: -1 to 1, 1 to 0
    Z2_number += ky_sign
  end
  println("[info] Bad point number is: $(bad_point_num)/$(ky_num)")
  return Z2_number % 2
end


### Write the data to files
function save_data_to_file(data, file)
  open(file, "w") do fwp
    writedlm(fwp, data)
  end
end


### Calcualte the Z2 number via wlison loop
function calc_Z2_via_wlison_loop(
  kx_beg, kx_end, kx_num,
  ky_beg, ky_end, ky_num, 
  occupy_num, data_file, data_type)
  ## Calculate the wilson loop
  k_path_y, wilson_loop_for_bands = 
    calc_wilson_loop(
      kx_beg, kx_end, kx_num, ky_beg, ky_end, ky_num, 
      occupy_num, data_file, data_type)
  ## Calcualte the Z2
  println("[do] Calcualting the Z2 number...")
  Z2_number = sum_to_Z2(ky_num, occupy_num, wilson_loop_for_bands)
  # Save Wilson loop data to files
  println("[do] Saving data to files")
  save_data_to_file(k_path_y, "KPATH.dat")
  save_data_to_file(wilson_loop_for_bands, "WLOOP.dat")
  return Z2_number
end


### TEST code. Remove follows if use this code as a lib.
Z2_number = calc_Z2_via_wlison_loop(-0.5, 0.5, 201, -0.5, 0.5, 51, 18, "w90.dat", "wannier90")
println(Z2_number)  