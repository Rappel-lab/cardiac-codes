#Three Dimensional Reconstruction to Visualize Atrial Fibrillation Activation Patterns on Curved Atrial Geometry

This folder contains the necessary sample data and code to generate the electrode shell, interpolate
the phase signal in 2D and 3D, identify phase singularities in 2D and 3D, and project the 3D
phase singularities back into a 2D mesh for comparison purposes. The folder 'Patient_Data' includes sample
EGM signals and NavX data. The 'DATA' folder includes additional information extracted from 'Patient_Data'
like the termination site or the position of the electrodes.

In order to successfully run the code, when using it in the Matlab console, the user should first add to path  
all folders and subfolders. Then, the described functionality can run the by calling the interp_method 
function as follows:

interp_method('ID-00', 1);

The function documentation includes additional information about the inputs and outputs. We hope that 
the reader is able to follow the different sections in this 'interp_method' function. This function
includes calls to the additional functions included in the folder. The reader can also check the
documentation of these functions if necessary.

###Raw data of several patients can be downloaded from the following google drive link: 
https://drive.google.com/drive/folders/1eYm_J4m0nVPZ48Ujhw3EUrcEoZH-Nhbb?usp=sharing
