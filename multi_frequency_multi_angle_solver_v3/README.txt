20241219, Zhaowei Dai
This is an updated version of v2 codes.

To-Do:
Add formatted documentations for newly added functions.
Change the indentation of functions to indent all functions.

Compared to last version multi_frequency_multi_angle_sovler_v2:
List of changes:
1. fom:
   coupling_efficiency_v2: calcualte the overlap integral inside the waveguide.
2. simulation: 
   add_absorber: The absorber regions are aligned with the simulation grid. The patch positions are on grid. Alignment is done for convergence test.
3. plot:
   deleted plotEonStructure_show_sim_settings and plotEonStructure, and replace them with a new function:
   plotFieldonStructrue_show_sim_settings_custom_colormap: use the variable to be plotted as the input, real(E), abs(E), abs(E).^2 instead of E
   and we can choose the colormap to be used, e.g. meep, and set the color limit use climit
   wavelengthToRGB: a function that returns the RBG value for a wavelength
   singleColorIntensityColormap: a function that generates the colormap for a single color with varying intensities
4. objective_function
   I put the mean_obejective in this folder, it's a very general function that do not rely on simulation settings
5. grating_cooupler: 
   compute_CE_and_gradients: this function is not general and is only useful for grating coupler optimization. But the structure of the code
   can be applied to other multi-angle multi-wavelength optimizations. You need to change the fom function for other foms.
   compute_CE_and_gradients_parfor: use parfor lopps instead of foor loops to parallelize frequencies
6. geometry:
   add_region: align with simulation grids, the 4 corners of the region are all on grids.
7. waveguide_mode:
   TE_symmetric_waveguide: Calculates TE propagation constant and field profiles for a symmetric waveguide 
   compute_wg_mode_and_power: Compute waveguide mode profile along a line monitor and power flow across the monitor for multiple wavelengths
   and multiple modes at each wavelength. 
8. source:
   GaussianBeam2D_EzHy: % Calcualate incident fields and power, provide Ez Hy function handles to evaluate fields on given position

   






   
   



