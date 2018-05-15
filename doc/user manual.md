# User Manual for hMC-GRT v1 (04/2018)

## A. ZEMAX Excitation 
1. Using the ZEMAX non-sequential mode, build your optical model for excitation path including your sample. 

2. Make an object of ‘Rectangular Volume’ mimicking a semi-infinite and homogenous medium of your sample and locate it as your sample position (typically z axis is an optical path direction in ZEMAX). Set up a proper refractive index of your sample. 

3. Place the ‘Ellipse’ underneath this sample surface (relative depth: - 1e-4 mm). This Ellipse will be your detector so make sure to create this large enough to capture all the incident rays. The captured rays by this ‘Ellipse’ go through refraction (incident angle change) and this new angle will be an actual incident angle of each photon launched into a sample for Monte Carlo simulation. 

4. To run the simulation, click Analysis > Ray Tracing > Racy Trace… (Ctrl +D) 

5. In the ‘Ray Trace Control ‘window, check ‘Save Rays’. Put the filename ‘N’-filename.dat (‘N’ is an object number of ‘Ellipse’ made at the step A.2. 

6. Make sure that # of analysis ray is 106. (This is the maximum number)  Run the simulation. Find the resulted ‘filename.dat’ in /Document/ZEMAX/Object/Sources/Sourcefile 

7. Run ‘SourceFileRead.exe’ to convert ASCII-type data file to text data file. For detailed information on how to use, please use a link below 

http://customers.zemax.com/os/resources/learn/knowledgebase/how-to-convert-a-binary-source-file-into-ascii 

8.  Once conversion completed, open the output text file in Microsoft Excel. Remove the first 4 rows containing file conversion information. Save this file (For explanation purpose, let’s say it is “zemax_input.dat”). Remember the number of rays for the following Monte Carlo simulation. This output text file contains seven columns of information per each ray hitting the ‘Ellipse’
col 1: x coordinate
col 2: y coordinate 
col 3: z coordinate 
col 4: 𝛼 cosine 
col 5: 𝛽 cosine 
col 6: 𝛾 cosine 
col 7: intensity

## B. Monte Carlo Simulation (For detailed information, please refer to an original user manual of MC simulation written by Karthik Vishwanath)

1.  Edit your “filename.tissue” file for your sample’s optical properties. Change the number of photons as the same with the number of rays in A.8

2.  Edit “trmc_2fls_globals.h” to modify the simulation options. Put the filename of the ZEMAX input (zemax_input.dat) and output array file (e.g. zemax_output.dat). Make sure the input array file name is the same as A.7

3.  Open a terminal program and run ‘Makefile’ to build an executable file named ‘trmc_mlmf’ 
```
/> make time 
```
  3.1. Uploaded “Makefile” was tested to compile and link the source codes using a MinGW complier (gcc, www.mingw.org) in Window 10 environment. 

  3.2. Linux environment (Ubuntu) was also tested. “Makefile” needs to be modified for a different complier option.

4.  Execute the simulation - 
```
 /> trmc_mlmf filename.tissue
```
5.  As a result of MC simulation, ‘zemax_ouput.dat’ file will be generated. This output text file contains eight columns of information per each photon exiting the surface of the sample 

col 1: x coordinate
col 2: y coordinate 
col 3: z coordinate 
col 4: 𝛼 cosine 
col 5: 𝛽 cosine 
col 6: 𝛾 cosine 
col 7: intensity
col 8: wavelength 


## C. ZEMAX Detection 

1. Open the “zemax_ouput.dat” in Microsoft Excel. 

2. We need to update this file to be compatible with the ZEMAX source file format. Insert the total number of photons in the first row and first column. Insert ‘2’ in the first row and second column. ‘2’ indicates ‘cm’ unit because MC simulation used ‘cm’ unit. 

3. Save and move this file in the folder (/Document/ZEMAX/Object/Source/Sourcefile)

4. In ZEMAX non-sequential model where your detection model is created,  

  -4.1 Insert “Source File” object, select “Data File” and open the file saved in C.2 

  -4.2 Change z-coordinate of this source file at the surface of your virtual sample, which mimics photons exiting your sample and will be a light source for ZEMAX simulation. 

  -4.3 Change the number of Analysis Rays to total number of photons (first row first column of excel sheet)

  -4.4 Change Power (Watts) to match Watts in File to make each ray 1W. 

  -4.5 Put the “Detector” object at the place where your actual detector (imaging sensor or detection fiber) is located. 

5. Run the ZEMAX simulation (Control D).  Select Save Rays. 

6. To see the results, go to Analysis -> Ray Tracing -> Ray Database Analysis -> Flux vs Wavelength
