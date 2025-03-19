
  Congratulation!  You now have a scanner calibration photo.
  Some important things you should keep in mind:

- This CD does contain the IT 8.7/2 reference file for the
  calibration target charge "R090104". Select the file "R090104.txt"
  as reference file when using the scanner profiling software.
  You can always download the latest reference files
  from http://www.targets.coloraid.de .

- Handle the targets carefully to avoid kink marks, scratches
  and fingerprints. Do not touch the surface.

- Remove the protection sleeve of the photo before scanning.
  Return the test target to it's protective enclosure immediately
  after use.

- Always protect the calibration photo from strong light.
  Store in a cool dry and dark place, preferable 70 F (21 C)
  and not exceeding 50% RH. Avoid sudden temperature changes
  as this can cause moisture on the photo surface.
  Use a dry lint free cloth to dry/clean the surface.
  DO NOT USE ALCOHOL or other chemicals to clean the surface
  of the photo.

- The colors of the target will change with time.  It is recommended
  to get a new target every 2 years if you want to avoid any
  quality loss. Please check out http://www.targets.coloraid.de
  for the latest information on your charge.

- You  can  use  the calibration photo with basicly any other
  proper  IT  8.7/2  software on PC/Mac/Amiga/Sun and other systems.
  So  far  the  target  has  been  tested  with ScanOpen (see "Extras"
  directory below!), ColorTune, CorelDraw, Silverfast, lcms and others.
  Note that the reference file is identical for all systems and
  thus can simply be transfered to other systems in order to get used.

- Users of the BasICColor profiling software having problems loading
  the reference file should copy the reference file into the
  reference files directory for Agfa targets of their BasICColor software.
  BasICColor may need additional an additional description file
  found in the Agfa directory to load the IT8 reference file.
  
- Check out http://www.coloraid.de for links to various scanner
  profiling software and other CMS tools.
  
- The "Extras" directory does contain additional files for your
  target:

  R090104W.txt :  Same as the normal reference file, but measurement
                  was done using a white instead of black backing.
                  The IT 8.7/2 standard requires the measurement of the
                  calibration target using a black backing.
                  Some (especialy consumer) scanners do have a white
                  backing. Use this special reference file if your
                  scanner does shipp with a white lid.

  R090104S.txt :  A smaller version of the normal reference file not
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

  Fault.txt    :  Some technical info on the measurement and
                  production of the target. See also the statistical
                  info saved in the normal reference file.

  R090104.hist :  This file is for use with the X-Rite ColorShop X program 
                  and does contain the spectral data of the target. Using 
                  the spectral data ColorShop is able to calculate a large 
                  number of color operations, display the gamut hull and 
                  more. Note that the file format only supports spectral 
                  data from 400 to 700nm. A fully functional but time 
                  limited demo version of the ColorShop X software for 
                  MacOS and Windows can be found on the X-Rite website: 
                  http://www.x-rite.com   For more information on the 
                  spectral data please read below.

- The spectral data files provided in the "Extras" directory are currently 
  mainly for use by experienced users and developers. Normal users 
  intending to profile their scanners do not need these files.

  The spectral data can mainly be used to calculate color data for different 
  observers or color spaces. For instance, if you want to know the color
  under D65 instead of D50 light, you can use ColorShop to calculate
  the values.

- Measurement is based on wavelengths from 380-780nm in 3nm intervals and 
  does meet the ISO 13655 standard. The 3nm data is interpolated to 10nm 
  data according to ISO 13655:1996(E) Annex A. For batch average measured 
  productions the spectral data available is the mean data of several 
  targets measured. XYZ tristimulus and other color spaces calculated from 
  the mean spectral data can differ from the mean color space values 
  listed in the IT 8.7 reference file. Status T density informations found
  in the reference file are according to ANSI CGATS.5-1993.

Wolf Faust - wfaust@coloraid.de
