# kinematic_corrections
A script containing the corrections presented in (Harborne et al., in prep 2019). These corrections can be used to remove the effects of beam smearing from observed kinematics - &#955;<sub>R</sub> and V/&#963;. 

Specify the input CSV file (-i) and the output CSV file (-o). The corrected kinematics will be added to the original file as column that reflect the observed kinematics (i.e. if "obs_lr" column is present, a "corr_lr" column will be added).

To run from the command line:
> python kinematic_corrections.py -i "path_to_your_input_CSV.txt" -o "path_to_your_output_CSV.txt"

This code accepts an input CSV file with columns for:
* Observed &#955;<sub>R</sub> ("obs_lr") <b>AND/OR</b> observed &#955;<sup>&#949;</sup><sub>R</sub> ("obs_elr") <b>AND/OR</b> observed V/&#963; ("obs_vsig").
* Observed ellipticity of the measurement region ("e").
* Sersic index ("n").
* Observed seeing conditions divided by the measurement radius i.e. PSF_fwhm/2.355)/Re  ("psf_over_re") <b>OR</b> the full-width-half-maximum of the PSF ("fwhm") <b>AND</b> the semi-major axis length radius of the  measurement region ("Re").
* The number of effective radii that the measurment radius encompasses i.e. Reff/Re ("Reff_fac").

The order of these columns is <b>not</b> important, but all must be present for the script to run without returning an error. 

If any column is missing, you will see an error or warning message. If "Reff_fac" is missing, the code assumes your kinematics have been measured within 1 effective radius, but will warn you of this in the output messages and produce an output file. If any other columns are missing, no output file will be produced. 

