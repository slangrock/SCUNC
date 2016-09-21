# SMELLIE Calibration Utensils for Normal scattering length Calculation - SCUNC

Framework to define the SMELLIE Monte Carlo simulation and carry out the Rayleigh scattering length of a water-filled SNO+ detector.

## Description

Light scattered in the detector arrives later than expected in any PMT. Therefore, the reconstruction and energy resolution of a detector are dependent on the knowledge of the scattering properties of the detector medium. The SMELLIE system was designed to measure the scattering properties of the SNO+ detector. The systematic uncertainty on the scattering length is dependent on how well the simulations of the system describe the data. Hence, this framework was developed to carry out both the scattering length measurement as well as the Monte Carlo definition of the system. The files contained within this framework are:

1. Monte Carlo description:
   * fibre_direction_measurement.py
   * get_angle_profiles.py
   * plot_angleprofiles.py
   * sys_uncert_angle.py
   * intensity_measurement.py
   * get_intensity.C

2. Scattering length measurement:
   * apply_cuts.py
   * cut_verification.py
   * plot_cuts.py
   * measure_scalingfactor.py

3. Additional scripts:
   * check_files.py
   * bias.py

4. Files containg classes and functions used by other scripts:
   * angle_measurement.py
   * data_corrections.py
   * define_histograms.py
   * fibre_handling.py 
   * plot_style.py

To run any of these files ROOT, the SNO+ RAT and a python environment have to be loaded. Please be aware that in the current state all scripts contain hardcoded directory paths which will have to be adjusted for different users. To run the scripts, two text files named fibre_positions_new.txt and cut_values.txt are needed, which are not included in this repository but can be made accessible upon request.

## Monte Carlo Definition

### Fibre direction

To measure the direction vector of a set of fibres run 

    $ python fibre_direction_measurement.py -f 'fibre_label' -o direction

The -f flag assigns which fibre is measured, the -o flag determines if a direction measurement should be carried out. The other -o option is 'systematics' and returns a set of random variables in which the uncertainties determined during the direction measurement can be applied.

### Angular Profile

To extract the angular profiles from a SMELLIE data run do

   $ python get_angle_profiles.py -m 'medium' -p 'n_pmts'

The script was written when only air-fill data was available. Thus, assigning 'air' to the -m flag returns a profile extracted straight from a data file. To carry out a water-fill scattering length measurement, a method was developed to convert the profiles from air data to a water-fill profile when 'water' is applied as -m flag. The profiles are saved into root files containing histograms with variable angular bins. The -p flag assignes how many PMT centers are covered by each angular bin. The bins of the profiles are defined using the AngleMeasurement class (angle_measurement.py). All measured profiles are also averaged and smoothed.

To create plots and convert the profiles into probability arrays as used in SMELLIE.ratdb do

   $ python plot_angleprofiles.py -f 'fibre_label' -w 'wavelength'

The profiles can be altered to evaluate systematic uncertainties arising from the profiles of the fibres. The altered profiles can be produced with

   $ python sys_uncert_angle.py -f 'factor' -t 'type'

Both the width and the tail contributions of the angular profiles need to be considered for systematic uncertainties. The -t flag allows to alter 'width' or 'tail' separately or 'both' at the same time. The -f flag attributes the factor by which the profile should be altered.

### Intensity Measurement

To measure the intensity of a SMELLIE run, the root files have to have the structure FIBRE_WAVELENGTH_INTENSITY.root The intensity of SMELLIE data can be measured by first running

    $ python intensity_measurement.py -d 'dir_name' -f 'fibre_label' -w 'wavelength' -t 'type'

Dependent if the 'type' is 'sim' or 'data' the function will save the fit parameters of a nhits vs intensity function or the nhits mean into a text file. The -d flag points to the directory of the root files. To use these values to achieve the intensity of the data in a format that can be used in the simulation, do

    $ root -l
    $ [0] .L get_intensity.C
    $ [1] main('fibre_label','wavelength')

This script saves the measured intensity and its uncertainty to a text file.

## Scattering Length Measurement

For the scattering length measurement, the SMELLIE root files have to be named FIBRE_WAVELENGTH_SCALINGFACTOR.root. To apply a cut selection to SMELLIE runs, run

   $ python apply_cuts.py -d 'dir_name' -f 'fibre_label' -w 'wavelength' -t 'type'

If the chosen type is 'sim' the cuts are applied to the photons in the MC branch of the root file. If the chosen type is 'data' the cuts are applied to the hits found in the EV branch of the root file. This script returns root files containing histograms defined by the DefineHistograms class in define_histograms.py and filled with the events found for each cut region, as well as a text file containing all events found in each cut region and the ratios of in-beam and scattered light.

To check how many scattered or noise photons are selected by each cut, run

   $ python cut_verification.py -d 'dir_name' -f 'fibre_label' -w 'wavelength' -t 'type'

with type being either 'tracks' or 'noise'.

The histograms filled by apply_cuts.py can be plotted using

   $ python plot_cuts.py -f 'fibre_label' -w 'wavelength' -r 'ratio'

with ratio being the scaling factor of the SMELLIE run.

The scaling factor is measured from the text files produced by apply_cuts.py using

   $ python measure_scalingfactor.py -d 'dir_name' -f 'fibre_label' -w 'wavelength'

with the -d flag pointing to the directory containing the text files.

## Additional Scripts

To check if the input root files in a directory are valid root files do

   $ python check_files.py -d 'dir_name'

The bias of the scaling factor measurement can be determined from the scaling factors returned by the measure_scalingfactor.py script using

   $ python bias.py -f 'fibre_label' -w 'wavelength' -t 'type'

The -t flag allows for 'fibre' or 'wavelength' specific measurements only, or for a bias measurement over all fibre-wavelength combinations ('full').