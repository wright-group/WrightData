# Measurement of Ultrafast Excitonic Dynamics of Few-Layer MoS<sub>2</sub> Using State-Selective Coherent Multidimensional Spectroscopy

Authors: Kyle J. Czech, Blaise J. Thompson, Schuyler Kain, Qi Ding, Melinda J. Shearer, Robert J. Hamers, Song Jin, and John C. Wright

This material is excerpted from a work that was published in ACS Nano, copyright © American Chemical Society after peer review. To access the final edited and published work see http://pubs.acs.org/doi/abs/10.1021/acsnano.5b05198 .

## Contents

1. [Movie](#movie) - the full three dimensional frequency frequency delay spectrum presented as the main result
2. [Absorbance](#absorbance) - absorbance spectrum of the MoS<sub>2</sub> thin film
3. [Reflective 2D Delay](#reflective-2d-delay) - two dimensional delay delay spectrum of the MoS<sub>2</sub> thin film in reflective geometry
4. [Transmissive 2D Delay](#transmissive-2d-delay) - two dimensional delay delay spectrum of the MoS<sub>2</sub> thin film in transmissive geometry
5. [Images](#images) - images of thin film

## Movie
The full three dimensional (w1, w2, d2) coherent multidimensional spectrum presented as the main result in this work.

### Location:
\\movie

### Data description:
The three dimensional data was collected as a series of w1 d2 wigner spectra. Each wigner spectra is recorded in its own .dat file. Note that a fairly extensive processing procedure was performed on the data prior to generating the plots presented in this work. This procedure is described in the published work. In addition, you may replicate the processing procedure using our data processing package WrightTools and the `workup.py` script located in the same folder as this README. WrightTools is free and open access and can be downloaded from [github](https://github.com/wright-group/WrightTools).

### File description:
The data is stored in a series of .dat files generated by COLORS (COntrol for Lots Of Research in Spectroscopy). Dat files are tab-delimited text files. Each row is an acquisition corresponding to a given set of hardware positions (a pixel in the experimental space). The columns are as follows. Columns containing data presented in this experiment are indicated in bold. In the internal Wright group labeling scheme, this is COLORS dat version 2.

1. acquisition index
2. **OPA1 color (nm) (w1)**
3. polarization of OPA1 (0 for vertical, 1 for horizontal)
4. **OPA2 color (nm) (w2)**
5. polarization of OPA1 (0 for vertical, 1 for horizontal)
6. OPA3 color (nm)
7. polarization of OPA3 (0 for vertical, 1 for horizontal)
8. monochromator color (nm)
9. array H index
10. array V color (nm)
11. reference delay position (fs)
12. reference delay zero position (mm)
13. delay 1 position (fs)
14. delay 1 zero position (mm)
15. **delay 2 position (fs) (d2)**
16. delay 2 zero position (mm)
17. **DAQ analog input 0 (V) (signal)**
18. DAQ analog input 1 (V)
19. DAQ analog input 2 (V)
20. DAQ analog input 3 (V)
21. DAQ analog input 4 (V)
22. array value (arbitrary units)
23. OPA motor 0 position (us)
24. OPA motor 1 position (us)
25. OPA motor 2 position (us)
26. OPA motor 3 position (us)
27. OPA motor 4 position (us)
28. OPA motor 5 position (us)
29. OPA motor 6 position (us)
30. calculation 0
31. calculation 1
32. calculation 2
33. calculation 3
34. calculation 4
35. calculation 5

Some columns were not used during the collection of this data. In particular, the OPA3, array, OPA motor, and calculation columns do not contain information that is relevant to this work in any way.

## Absorbance

TO DO

## Reflective 2D Delay

TO DO

## Transmissive 2D Delay

TO DO

## Images

Images taken of the MoS<sub>2</sub> thin film to show homogeneity and smoothness.

### Location:
\\images of thin film

### File description:
.png