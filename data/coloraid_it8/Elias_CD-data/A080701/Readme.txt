  Congratulation!  You now have a scanner calibration target.
  Some important things you should keep in mind:

- The ANSI IT 8.7/1 target is for use with the Agfa RSX, RSX II
  and Ct Precisa films. Using other films usualy does result in
  visible color faults.

- This floppy or CD does contain the IT 8.7/1 reference file(s) for the
  calibration target charge "A080701". Make sure that the charge printed on 
  your target slide frame does match the reference file. For instance, 
  choose the file "A080701.txt" in your profiling software if the charge 
  printed on the slide frame is "A080701".

- Handle the targets carefully to avoid kink marks, scratches
  and fingerprints. Do not touch the surface.

- Always protect the calibration slide from strong light.
  It is strongly recommended to store the target and your slide
  films in a cool, dry and dark place.

  A good and more detailed guidance in handling and storage of
  photographic materials can be found in Kodak's Publication
  No. E-30, "Storage and Care of Kodak Photographic Materials - 
  Before and After Processing". The document is available on
  Kodaks website http://www.kodak.com.

- The colors of the target will change with time.  It is recommended
  to get a new target after 2 years if you want to avoid any
  quality loss. Please check out http://www.targets.coloraid.de
  for the latest information on your charge.

- You  can  use  the calibration target with basicly any other
  proper  IT  8.7/1  software on PC/Mac/Amiga/Sun and other systems.
  Check out http://www.coloraid.de for links to various free scanner
  profiling software and other CMS tools.

- Users of the BasICColor profiling software having problems loading
  the reference file should copy the reference file into the
  reference files directory for Agfa targets of their BasICColor software.
  BasICColor may need additional an additional description file
  found in the Agfa directory to load the IT8 reference file.

- The color gamut of film is usualy larger than most output hardware
  can handle. Also not all scanner profiling programs can deal with large 
  color gamuts of slide films equaly well. Highly saturated colors possible 
  with film slides are sometimes not scanned or reproduced correctly. If 
  you run into huge color faults with highly saturated colors, try to 
  isolate the faulty color profile involved. It might be the output and not 
  the scanner input profile used. In case the scanner profile is faulty, 
  try a different profiling software that has been well tested with 
  slide films.

- The "Extras" directory does contain additional files for your
  target:

  A080701S.txt :  A smaller version of the normal reference file not
                  containing any statistical information on the target
                  production, resulting in a smaller file.
                  Use this reference file if your profiling
                  software does not accept the standard reference file.
                  Currently there seem to by only two profiling
                  programs not able to accept the standard reference
                  files. A fault in some versions of the
                  Heidelberg CPS ScanOpen software do cause a crash
                  when reading the large reference file.
                  Agfa's ColorTune users should also use this smaller
                  file in case the normal file is not accepted
                  by the profiler.

  A080701.eps  :  This file is for use with the X-Rite ColorShop V2 program
                  and does contain the spectral data of the target. Using 
                  the spectral data ColorShop is able to calculate a large 
                  number of color operations. Note that the file format 
                  only supports spectral data from 390 to 690nm. Tests 
                  indicate that the error introduced by this limitation is 
                  small (roughly dE 0.2).  Please read below for more
                  information about the spectral data.

  A080701.xls :   A tab delimited file that loads nicely into the Excel
                  spreadsheet program. The file contains the spectral data
                  of the targets and various color values calculated from
                  the spectral data. Please read below for more information
                  about the spectral data.

  A080701.cgt :   This file does contain the spectral data of the targets
                  using the ANSI CGATS.5-1993 format. Despite being very
                  similar to the IT 8.7 reference file, you should not use 
                  this file for generating scanner profiles unless noted 
                  differently. Some of the color data stored in this file 
                  may differ from the normal IT 8.7 reference file. Please 
                  read below for more information about the spectral data.

  A080701.cxf :   This CxF formated file does contain the spectral data of 
                  the targets.  The XML  based CxF  format is  mainly used 
                  by newer Gretag -     Macbeth     soft-    and hardware 
                  products. Eye-One Share (http://www.i1color.com/freeware),
                  ProfileMakers MeasureTool and ColorLab (both available 
                  from http://www.gretagmacbeth.com) are programs available 
                  online at  least  as  partly  functional  demo version 
                  and  do allow loading/displaying the  file.
                  All CxF compatible programs mentioned seem to 
                  have limitations and do not make full  use of the data 
                  similar to ColorShop. The data is limited to 380-730nm.

- The spectral data files provided in the "Extras" directory are currently 
  mainly for use by experienced users and developers. Normal users 
  intending to profile their scanners do not need these files.

  The files can mainly be used to calculate color data for different 
  observers or color spaces. For instance, if you are used to work with 
  density data than you can use X-Rites ColorShop to display the density 
  values for each patch.

- Measurement is based on wavelengths from 380-780nm in 3nm intervals and 
  does meet the ISO 13655 standard. The 3nm data is interpolated to 10nm 
  data according to ISO 13655:1996(E) Annex A. For batch average measured 
  productions the spectral data available is the mean data of several 
  targets measured. XYZ tristimulus and other color spaces calculated from 
  the mean spectral data can differ from the mean color space values 
  listed in the IT 8.7 reference file. Status A density informations found
  in the reference file are according to ANSI CGATS.5-1993.

Wolf Faust - wfaust@coloraid.de
