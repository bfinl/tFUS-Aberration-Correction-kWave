# kPR: A k-Wave based Phase-Reversal Method for tFUS Skull Aberration Correction

Transcranial focused ultrasound (tFUS) offers significant promise for precise, non-invasive neuromodulation, with potential applications in treating neurological disorders and advancing neuroscience research. However, skull-induced acoustic aberrations critically hinder the effectiveness of tFUS by distorting focal precision and reducing the energy delivered to targeted brain regions.

In this study, we introduce and validate a novel phase-reversal aberration correction method (kPR) designed specifically for a 128-element phased-array ultrasound transducer. Our approach employs individual head models derived from magnetic resonance (MR) data of 22 subjects and sophisticated k-Wave simulations to address these aberrations. Results show substantial improvements, including a 98.70% increase in overlap volume between the ultrasound focus and targeted brain region, a 14.36% reduction in axial positioning errors, a 21.53% decrease in focal peak positioning errors, and a 17.58% enhancement in targeted energy delivery compared to cases without correction. In addition, we demonstrate the superiority of our approach over traditional ray-based correction methods through comprehensive comparative analyses.

This advancement significantly elevates the potential of tFUS for precise, personalized neuromodulation therapies, paving the way for safer, more effective clinical interventions and robust investigative neuroscience studies.

The codes we share here are developed and evaluated in the following study, and is provided without warranty. Users should use at their own risks in accordance with the [licensing requirements](#license).

> Z. Li, K. Yu, J. Kosnoff, and B. He, "Improving Targeting Specificity of Transcranial Focused Ultrasound in Humans Using a Random Array Transducer: A k-Wave Simulation Study," *bioRxiv*, 2025. DOI: 10.1101/2025.04.25.650630

If you use our kPR method or any part of our codes in your own work, please acknowledge us by citing the above manuscript.

## Platform and Dependencies

Tested on Linux (Ubuntu 20.04.3 LTS, MATLAB R2021a) and on Windows 11 (MATLAB R2023a and R2024b). Works with both NVIDIA GPU and CPU-only platforms:

- For CPU-only k-Wave simulations, you only need to use k-Wave ver 1.4 or higher ([download](http://www.k-wave.org/download.php)) and select `model = 1` under the DEFINE LITERALS sections in:
  - `./H275_3D_human_phase_correction.mlx`
  -  `./H275_3D_human_phase_estimation_local.m` 
  -  `./generate_source.mlx`

- For GPU-accelerated simulations, you can consider using our compiled C++/CUDA binaries provided under `./k-Wave binaries/` (by replacing the original binaries provided in the k-Wave toolbox under `k-Wave/binaries/`) and select `model = 4`, which allows you to use single NVIDIA GPU up to the *Ada Lovelace* architecture (RTX 40 series). In addition, you have the option of using the MATLAB GPU code by selecting `model = 2`.
- *Note: we do not recommend using `model = 3` (compiled C++ code, CPU-only) since it is not widely tested by us.*

## Instructions

1. Prepare subject head CT (or pseudo-CT) and MR data (preferably in `nifit` format otherwise you need to do conversion first). If you only have MR data of subjects, we recommend you follow the instructions under `./MR to pCT/` to convert it to pseudo-CT (pCT).
2. Follow the instructions under `./Brain MR Segmentation/` to obtain the segmented region of interest (V5L in our work). This would be used for quantitative evaluations and statistical calculations later on.

3. Use the `./skull_brain_visualization.mlx` script to visualize the skull (CT or pCT) and brain (MR) in a 3D space. Also, you can further determine the thresholding value for skull/brain binary mask generation using this script.

4. (Optional) Use the `./generate_source.mlx` script to create time varying source (CW signal in our example code) for single element and save them into separate files. This step can help you greatly reduce the total computational time (by hours) when using multi-element transducers (e.g., 128 or 256 elements) since the phase delay estimation will be done on a element-by-element basis. *Note: You should set parameters properly in each section before running the script.*
5. Set parameters properly in `./H275_3D_human_phase_estimation_local.m` and run the script to estimate phase delay values for each element of your transducer. This script will automatically save not only the estimated phase delay values, but also save the raw sensor data and figures of the waveforms at the target location. You may modify the code to turn off such extra functions if necessary. *Note: By default, this piece of script utilizes the pre-generated source matrices mentioned in the previous step. You may need to modify the code if you do not want to generate source beforehand.*
6. When phase delay estimation is done, use the `./H275_3D_human_phase_correction.mlx` script (be sure to set parameters properly also) to get your final results by firing all elements with the estimated phase values applied to the source signal. You may run the code with or without applying the estimated phase values for comparison purpose.
7. (Optional) Use the `./data_analysis.mlx` script for significance analysis (e.g., t-test) and generate box plots for multiple subjects' data you obtained.

## Acknowledgements

This work was supported in part by National Institutes of Health grants R01NS124564, U18EB029354, T32EB029365, RF1NS131069, R01NS096761, and R01NS127849-01A1, awarded to Dr. Bin He, Carnegie Mellon University.

## License

- **This project’s own code** is released under the **GNU Lesser General Public License (LGPL) v3**.  
  See [LICENSE](./LICENSE).
- **Bundled third-party code/compiled binaries** in:
  - `./k-Wave binaries/` — also under **LGPL v3**.

*Note: All modifications you make to the LGPL’d parts stay under LGPL v3.*