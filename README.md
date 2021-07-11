# Watta – matlab code for detecting melt water lakes![image](https://user-images.githubusercontent.com/47567083/125211734-b15a6e00-e265-11eb-9048-7300971fcb65.png)

Watta – a matlab code for detecting melt water lakes

This folder contains the code used in Datta, R. T. and Wouters, B.: Supraglacial lake bathymetry automatically derived from ICESat-2 constraining lake depth estimates from multi-source satellite imagery, The Cryosphere Discuss. https://doi.org/10.5194/tc-2021-4, in review, 2021.

The code is written in Matlab and uses the ATL03 and ATL03 file of the track to be processed as input. The code is initiated by running ‘run_ATL03_surfaceDet_Watta.m’.  Filenames and locations of the input and output directories should be defined in this file. The latitude extents of the part of the track to be processed should also be defined here. Parts of the code can be run in parallel, as indicated in the comments in the code (run_ATL03_surfaceDet_Watta.m and post_reset_kerndensity_lakes.m) 

We are currently working on moving the code to an open-source programming language. We therefore don’t foresee any major changes to the Matlab code, which is provided as is.

Pleases see pdf documents for a flow diagram explaining both steps, beginning with the process of developing the initial depth algorithm (Watta_code_DepthEst.pdf) and subsequent smoothing and feature extraction (Watta_code_SurfaceChar.pdf).

![image](https://user-images.githubusercontent.com/47567083/125211749-c0412080-e265-11eb-9080-6ccd43374b59.png)
