# BUDA-EPI

Demonstration code for the MRM paper, ***Liao, C., et al. "Distortion-free, high-isotropic-resolution diffusion MRI with gSlider BUDA-EPI and multicoil dynamic B0 shimming." Magnetic Resonance in Medicine (2021).

The datasets presented in this repository are from a BUDA-EPI acquisition with blip up and down shots. 

Written by Congyu Liao. Please feel free to contact me (cyliao@stanford.edu) if there's any question I can help answer.

**Data download:
```
The download link of BUDA-EPI datasets is: https://drive.google.com/file/d/1DAeTW2iLgQx8tjZT3QqZ7F3vvt6Vwt2Q/view?usp=sharing
```

The data size is ~ 20 GB, and the scan protocol is similar to the protocol of human connectom project (HCP):

```
FOV: 220 mm×220mm, matrix size: 176×176
In-plane resolution: 1.25 mm ×1.25mm
Slice thickness: 2mm; 
Number of slices: 57
Number of coils: 32
TR/TE = 2800/77 ms
Effective echo-spacing = 0.37 ms
Multi-band factor: 3 with CAPI-shift FOV/2 per shot
In-plane acceleration factor :2 per shot
partial Fourier factor: 0.75
BUDA acquisition with AP shot kyshift 1 and PA shot kyshift 0
Two shells: b= 1000 with 32 diffusion directions and b=2500 with 64 diffusion directions
```
**Tested with:

    MATLAB 2015b
    BART 0.5.00 (https://zenodo.org/record/3376744)
    FSL 6.0.4 (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

**Requirements:

    This repository requires BART and FSL to be installed. More details can be found in "instructions.pdf" in the folder
    
    



