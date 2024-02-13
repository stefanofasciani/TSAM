
=======================================
TSAM - Timbre Space Analyzer and Mapper
=======================================

==============================================================================

The Timbre Space Analyzer & Mapper (TSAM) is a collection of MAX patches 
and MATLAB functions for the analysis, modeling and malling of the
timbre of sound synthesizers.
The TSAM can be obtained at http://stefanofasciani.com/tsam.html
TSAM Copyright (C) 2016 Stefano Fasciani, University of Wollongong
Inquiries: stefanofasciani@stefanofasciani.com

==============================================================================


============
REQUIREMENTS
============

The TSAM is implemented in Cycling '74 Max 7 and MATLAB (OSX only).

The TSAM has been developed and tested on the following platform:
OSX 10.7.5
MAX 7.2.4 (32-bit)
MATLAB 7.12.0 (R2011a)


The following Cycling '74 Max 7 libraries are required:
(download the externals and include these in your Max path)

- FTM http://ftm.ircam.fr/index.php/Main_Page
- SFA-MaxLib http://stefanofasciani.com/?p=865
- CNMAT Max/MSP Externals (Opensoundcontrol and OSC-route only)http://cnmat.berkeley.edu/downloads
- Max Sound Box (ircamdescriptor~.mxo only) http://forumnet.ircam.fr/product/max-sound-box-en/


The following MATLAB libraries/toolboxes are required :
(download the files.m and include these in your MATLAB path)

- Instrument Control Toolbox http://www.mathworks.com/products/instrument/
- Statistics and Machine Learning Toolbox http://www.mathworks.com/products/statistics/
- Signal Processing Toolbox ttp://www.mathworks.com/products/signal/

============
INSTRUCTIONS
============

The TSAM.maxpat and TSAM_3DOGL.maxpat are equivalent. The only difference is in the timbre space visualization.
The TSAM.maxpat provides a 2D timbre space visualization, which is less CPU demanding.
The TSAM_3DOGL.maxpat provides a 3D timbre space visualization (OpenGL), which is more CPU demanding.
Both TSAM.maxpat and TSAM_3DOGL.maxpat provides 2D and 3D mapping capabilities.

Run the TSAM.maxpat (or TSAM_3DOGL.maxpat) and run the TSAM.m script in MATLAB.
The TSAM.m script communicates with the TSAM.maxpat via OSC.
To terminate the TSAM.m push the "Stop Engine" button in the TSAM.maxpat main tab.
When TSAM.m is running "Engine Stopped" changes state to "Engine Running".
The TSAM.maxpat request to the TSAM.m sctipt the mapping computation and scoe computation
(only when pressing the buttons "Compute Mapping" and "Compute Score" in the TSAM.maxpat main tab.
For any other operation the TSAM.maxpat do not require the TSAM.m to run in background (e.g. for analysis 
or real time mapping).

In MATLAB you can compile the TSAM.m into a matlab executable and run it from the terminal.
To compile in MATLAB use the command below.

mcc -m -v -I full_path_to_TSAM_MATLAB_folder -o compiled full_path_to_TSAM_MATLAB_folder/TSAM.m

Matlab will generate automatically a shell script run_compiled.sh, then you can run the TSAM.m runing 
sudo run_compiled.sh from the terminal (you need superuser credentials)

