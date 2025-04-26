# kPR: A k-Wave based Phase-Reversal Methods for Skull Aberration Correction

kPR is a wave theory-informed k-Wave based phase-reversal method to tackle skull-induced aberrations for transcranial focused ultrasound (tFUS) applications. Follow the instructions below to estimate and apply phase delay values to a multi-element array transducer in k-Wave simulations.

Be advised, the codes are provided as a service to the scientific community and are to be used at the user's own risk in accordance with the [licensing requirements](#license).

## Platform and Dependencies

Tested on Linux (Ubuntu 20.04.3 LTS, MATLAB R2021a) and on Windows 11 (MATLAB R2023a and R2024b). Works with both NVIDIA GPU and CPU-only platforms:

- For CPU-only k-Wave simulations, you only need to use k-Wave ver 1.4 or higher ([download](http://www.k-wave.org/download.php)) and select `model = 1` under the DEFINE LITERALS sections in:
  - `./H275_3D_human_phase_correction.mlx`
  -  `./H275_3D_human_phase_estimation_local.m` 
  -  `./generate_source.mlx`

- For GPU-accelerated simulations, you can consider using our compiled C++/CUDA binaries provided under `./k-Wave binaries/` (by replacing the original binaries provided in the k-Wave toolbox under `k-Wave/binaries/`) and select `model = 4`, which allows you to use single NVIDIA GPU up to the *Ada Lovelace* architecture (RTX 4000 series). In addition, you have the option of using the MATLAB GPU code by selecting `model = 2`.
- *Note: we do not recommend using `model = 3` (compiled C++ code, CPU-only) since it is not widely tested by us.*

## Instructions

1. Prepare subject head CT (or pCT) and MR data (preferably in `nifit` format otherwise you need to do conversion first). If you only have MR data of subjects, we recommend you follow the instructions under `./MR to pCT/` to convert it to pCT.
2.  Follow the instructions under `./Brain MR Segmentation/` to obtain the segmented region of interest (V5L in our work). This would be used for quantitative evaluations and statistical calculations later on.

3. Use the `./skull_brain_visualization.mlx` script to visualize the skull (CT or pCT) and brain (MR) in a 3D space. Also, you can further determine the thresholding value for skull/brain binary mask generation using this script.

4. (Optional) Use the `./generate_source.mlx` script to create time varying source (CW signal in our example code) for single element and save them into separate files. This step can help you greatly reduce the total computational time (by hours) when using multi-element transducers (e.g., 128 or 256 elements) since the phase delay estimation will be done on a element-by-element basis. *Note: You should set parameters properly in each section before running the script.*
5. Set parameters properly in `./H275_3D_human_phase_estimation_local.m` and run the script to estimate phase delay values for each element of your transducer. This script will automatically save not only the estimated phase delay values, but also save the raw sensor data and figures of the waveforms at the target location. You may modify the code to turn off such extra functions if necessary. *Note: By default, this piece of script utilizes the pre-generated source matrices mentioned in the previous step. You may need to modify the code if you do not want to generate source beforehand.*
6. When phase delay estimation is done, use the `./H275_3D_human_phase_correction.mlx` script (be sure to set parameters properly also) to get your final results by firing all elements with the estimated phase values applied to the source signal. You may run the code with or without applying the estimated phase values for comparison purpose.
7. (Optional) Use the `./data_analysis.mlx` script for significance analysis (e.g., t-test) and generate box plots for multiple subjects' data you obtained.

## Cite This Work

If you use our kPR method or any part of our codes in your own work, please acknowledge us by citing the following paper and the repository.

> Li Z, Yu K, Kosnoff J, He B: "Improving Targeting Specificity of Transcranial Focused Ultrasound in Humans Using a Random Array Transducer: A k-Wave Simulation Study." (2025).

Please also consider citing k-Wave (see their [website](https://www.k-wave.org/) and [GitHub repository](https://github.com/ucl-bug/k-wave?tab=readme-ov-file) for details).

## Acknowledgements

This work was supported in part by National Institutes of Health grants R01NS124564, U18EB029354, T32EB029365, RF1NS131069, R01NS096761, and R01NS127849-01A1, awarded to Dr. Bin He, Carnegie Mellon University.

## License

- **This project’s own code** is released under the **GNU Lesser General Public License (LGPL) v3**.  
  See [LICENSE](./LICENSE).
- **Bundled third-party code/compiled binaries** in:
  - `./k-Wave binaries/` — also under **LGPL v3**.

*Note: All modifications you make to the LGPL’d parts stay under LGPL v3.*