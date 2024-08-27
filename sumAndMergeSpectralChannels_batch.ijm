#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".nd2") suffix

// Input: 8-channel spectral image with transmitted in the last channel
// Output: 32-bit 3-channel image:
//		Channel 1 = sum of input channels 1-3
//		Channel 2 = sum of input channels 5-7 (channel 4 of input is unused)
//		Channel 3 = input channel 8, converted to 32 bit for merging
// Theresa Swayne, Columbia University, 2024 for Hapshepsut Jackson and Liz Chang 

setBatchMode(true);
run("Bio-Formats Macro Extensions");

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {

	print("Saving to: " + output);
	print("Processing: " + input + File.separator + file);
	path = input + File.separator + file;
    run("Bio-Formats", "open=&path");
	
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
	
	// create a new 3-channel image
	run("Merge Channels...", "c3=Short c2=Long c4=Trans create");
	selectImage("Composite");
	saveName = basename+"_merge.tif";
	
	saveAs("tiff", output + File.separator + saveName);
	close(); 


}
