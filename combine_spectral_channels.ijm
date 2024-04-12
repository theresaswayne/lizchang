// ImageJ Macro Script to merge spectral channels
// Input: 8-channel spectral image with transmitted in the last channel
// Output: 32-bit 3-channel image:
//		Channel 1 = sum of input channels 1-3
//		Channel 2 = sum of input channels 5-7 (channel 4 of input is unused)
//		Channel 3 = input channel 8, converted to 32 bit for merging
// Theresa Swayne, Columbia University, 2024 for Hapshepsut Jackson and Liz Chang 
// TO USE: Open a spectral image. Run the script. 
// 		You can immediately proceed to running the GP calculation, or save the resulting image for later use.

// ---- Get image information ----
id = getImageID();
title = getTitle();
dotIndex = indexOf(title, ".");
basename = substring(title, 0, dotIndex);
extension = substring(title, dotIndex);
getDimensions(width, height, channels, slices, frames);

// ---- Process image ----

// Project short wavelength channels 1-3
run("Z Project...", "stop=3 projection=[Sum Slices]");
selectImage("SUM_"+title);
rename("Short");

// project long wavelength channels 5-7
selectImage(title);
run("Z Project...", "start=5 stop=7 projection=[Sum Slices]");
selectImage("SUM_"+title);
rename("Long");

// pull out transmitted image channel 8
selectImage(title);
run("Duplicate...", "duplicate range=8");
selectImage(basename+"-1"+extension);
run("32-bit");
rename("Trans");

// create a new 3-channel image for input
run("Merge Channels...", "c1=Short c2=Long c4=Trans create");
selectImage("Composite");
rename(basename+"_merge.tif");

