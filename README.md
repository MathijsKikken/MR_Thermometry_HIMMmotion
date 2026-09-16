# MR thermometry of RF heating in the human brain at 7T using the harmonical initialized multi-echo model (HIMM) with SVD-based motion correction

***
#### M.W.I. Kikken, B.R. Steensma, E.F. Meliadò, C.A.T. van den Berg, and A.J.E. Raaijmakers
#### Computational Imaging Group for MRI Therapy & Diagnostics, University Medical Center Utrecht
***

### Overview
This code repository implements an algorithm script to separate changes in phase into temperature, B0 field drift, head motion, and respiratory motion. 
Unfortunately, head motion results in additional fluctuations in phase, which will therefore also affect the proton resonance frequency shift signal.
Luckily, the change in the field due to motion is linear with the type and magnitude of motion. If we therefore know (1) how the head moves over time (during dynamic scanning) 
and (2) how the field changes with every type of motion (rotations around the AP-axis, rotations around the LR-axis, rotations around the FH-axis, and FH translations), we could correct for motion-induced field fluctuations.

An article detailing the work has been submitted for publication to a peer-reviewed journal. Use of this algorithm is allowed given proper citation of the original work and in accordance with the license file.

***

### Data availability
Note that the code requires dynamic data of either a baseline (at 0% SAR) or a heating (at 100% SAR) scan. To correct for motion, the motion parameters and characteristic motion-induced field maps are also necessary.
Data is acquired in 3 orthogonal slices, from which the motion parameters are estimated. 
Additionally, some extra scans are performed in which the subjects deliberately move their head. Characteristic motion-induced field maps are extracted from these scans using singular value decomposition. Code to reconstruct these motion-induced maps from deliberate motion data can be found in the subfolder "GetMotionMaps".
Please contact the corresponding author to get access to an example dataset of dynamic baseline, heating, and motion scans. Corresponding measured rotation/translation parameters will also be included.

***

### Previous algorithms
The algorithm provided in this repository is based on the Multi-echo MR thermometry algorithm using iterative separation of baseline water and fat images with l1 regularization. 
This method was proposed by Poorman et al. (https://github.com/poormanme/waterFatSeparated_MRThermometry) and published in Magnetic Resonance in Medicine [1].

Our group focuses on the measurements of RF heating. RF heating distributions are more diffusive, so the principle of l1 regularization does not work anymore 
(as that requires a localized heating hotspot, with very different spatial distribution from the drift field). 
The l1 regularization was removed from the code. In order to still separate the drift fields from the temperature increases, we proposed to initialize the drift fields using near-harmonic 2D reconstruction (proposed by Salomir et al. [2]). 
Near-harmonic 2D reconstruction used the fat layer (not susceptible to temperature-induced phase changes) that surrounds the anatomy to estimate the drift field in the whole region of interest. 
This drift field is decomposed into spherical harmonics and the corresponding coefficients are iteratively optimized together with the spatial temperature map and the coefficients of the library.

***

### Need for motion correction
This worked well for stationary anatomical regions such as limbs [3]. However, the head is subject to motion: besides B0 field drift and temperature, the phase over time will therefore also change due to cardiac, respiratory, and head motion. 
The previously proposed HIMM MRT method is therefore extended with a motion-compensation scheme (to compensate for head motion) and a library of reference images (to compensate for respiratory motion).

- First, the algorithm determines which of the reference dynamics in the library is most similar to the dynamic image being evaluated.

- Then, the motion-induced B0 field is estimated based on the measured motion parameters and the field disturbances originating from head motion. These field disturbances are fitted from pre-scans in which volunteers are instructed to deliberately move their head in specified directions. More details on this fit are provided in the folder "GetMotionFields".

- Subsequently, an initial estimate of the B0 field drift is acquired through near-harmonic 2D reconstruction. 
The skull is insensitive to temperature-induced phase changes and can therefore be used to provide a good initial estimate of the B0 field drift in the transverse plane (because the skull completely surrounds the brain in this plane).

- After this initial estimate of the B0 field drift was acquired. 
The algorithm jointly updates the reference dynamics, the motion-induced field, and the temperature map to make sure the model best resembles the dynamic currently being evaluated. 

This algorithm code uses a Water/Fat separation to get the water, fat, static off-resonance, and transverse relaxation rate components of the reference dynamics. 
These components are extracted from multi-echo data through the water/fat separation algorithm using graph cuts (Hernando et al. [4]). This algorithm can be accessed through the ISMRM fat/water toolbox (https://www.ismrm.org/workshops/FatWater12/data.htm).

***

The following files are relevant for utilization of the algorithm:
* __RunHybridHIMMSVD.m__: Script to run algorithm
* __OptimizeHIMMSVD.m__: Main function script
* __x_update.m__: Sub-function scripts

***

[1] Poorman, M. E., Braškutė, I., Bartels, L. W., & Grissom, W. A. (2019). Multi‐echo MR thermometry using iterative separation of baseline water and fat images. Magnetic resonance in medicine, 81(4), 2385-2398.

[2] Salomir, R., Viallon, M., Kickhefel, A., Roland, J., Morel, D. R., Petrusca, L., ... & Gross, P. (2011). Reference-free PRFS MR-thermometry using near-harmonic 2-D reconstruction of the background phase. IEEE transactions on medical imaging, 31(2), 287-301.

[3] Kikken, Mathijs WI, et al. "Multi‐echo MR thermometry in the upper leg at 7 T using near‐harmonic 2D reconstruction for initialization." Magnetic Resonance in Medicine 89.6 (2023): 2347-2360.

[4] Hernando, Diego, et al. "Robust water/fat separation in the presence of large field inhomogeneities using a graph cut algorithm." Magnetic Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine 63.1 (2010): 79-90.
