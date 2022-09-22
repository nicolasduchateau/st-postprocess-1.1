Code for the post-processing of speckle-tracking data, available under the license [CeCILL-B](http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html)

# Speckle-tracking post-processing v1.1 - Matlab

Release 1.1 = MATLAB scripts
(slight bugs corrected from version 1.0)

Author: Nicolas Duchateau (nicolas.duchateau@creatis.insa-lyon.fr)

Date: January, 2014

Links to the corresponding publications at: <br/> https://www.creatis.insa-lyon.fr/~duchateau/#publications

------------------------------------------------------------------------------------------------------------------------
**NOTICE:**

This code is made open-access. Comments and bug reports are welcome, as well as feedback on its possible improvements.

**Published reports of research using this code (or a modified version) should cite the following article that describes the method:**

*Duchateau N, De Craene M, Piella G, et al. A spatiotemporal statistical atlas of motion for the quantification of abnormalities in myocardial tissue velocities. Medical Image Analysis, 2011;15(3):316-28.*
https://doi.org/10.1016/j.media.2010.12.006

**The present MATLAB implementation is the one detailed in:**

*Duchateau N, De Craene M, Pennec X, et al. Which reorientation for the atlas-based comparison of motion from cardiac image sequences? In: Proceedings of Spatio-Temporal Image Analysis for Longitudinal and Time-Series Image Data, MICCAI'12 Workshop. Springer LNCS, 2012;7570:25-37.*
https://doi.org/10.1007/978-3-642-33555-6_3

------------------------------------------------------------------------------------------------------------------------
**IMPORTANT NOTE:** The data reading part is designed for data exported from ECHOPAC (GE Healthcare, Milwaukee, WI), using the "store full trace option". The user should adapt this part to the data format exported from other software.

**ARCHIVE CONTENT:**

- "Code" folder = MATLAB scripts. Scripts should be launched in the following order:
1) a1_ReadExportedData.m
2) a2_SpatioTemporalResampling.m
3) a3_PlotOutput.m

Functions "get_extension.m", "getECGnormalizedScale.m" and "getSTdataXY.m" are annex functions automatically called in the previous scripts, but do not need to be launched.

- "Data" folder = Sample data to test the scripts. Should be updated by the user:
1) "DEMO_DATA.xls" = Excel file containing the list of subjects to be processed + the information needed by the processing. See details further on.
2) Folders = one folder per subject. Each folder may contain several ".CSV" files, each of them corresponding to a single sequence.

------------------------------------------------------------------------------------------------------------------------
**RECOMMENDATIONS:**

A/ Data exportation: the "store full trace" option in ECHOPAC generates a .CSV at a given location of your disk ("Export" folder of the ECHOPAC directory). This .CSV file should be renamed according to the data stored in the .XLS subjects list and copied to the folders to be processed (see template from DEMO files).

B/ Convention for doing speckle-tracking: with this exporting option, data ordering within the .CSV file is sensitive to the way you do the speckle-tracking segmentations: at which location of the wall you do start and if you go clockwise or the opposite. Recommended protocol would be to start at basal-septal level (for 4CH views) or anteroseptal level (for SAX views) and first cover the septum.

C/ Filling the .XLS file. It should contain the following:

C.1/ file names (e.g. VOL / PIG / OFF / POST etc.)

C.2/ ECG events (mandatory if you want to temporally align data) = frame numbers read from the speckle-tracking interface. The present version expects: beginning of the cycle (onset of QRS), mitral valve closure, aortic valve opening, aortic valve closure, mitral valve opening, end of the cycle (onset of QRS).

C.3/ whether the data is 4CH or not, and if you need to flip it or not (normally you should not).

------------------------------------------------------------------------------------------------------------------------
**CONTENT OF THE MATLAB FILES:**

1) **a1_ReadExportedData.m** = reading data + getting displacement / velocity / strain / strain in Cartesian coordinates and in radial/longitudinal (or radial/circumferential if you're in short-axis) coordinates.

2) **a2_SpatioTemporalResampling.m** = until now, data has different number of points + instants for each processed subject. This resamples spatiotemporally, so that every subject has the same number of points + instants. This step also aligns events to a common temporal scale (you can specify which subject). Data is then ready for inter-subject comparison.

3) **a3_PlotOutput.m** = example for plotting data at a given location along the cycle. Vertical bars indicate the location of the involved ECG events (cf. Recommendation C.2)
