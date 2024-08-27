/* ImageJ macro for GP images analysis
   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms
   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus
   Nature Protocols 2011
   version July 2011
*/

/* This macro calculates GP images from Laurdan and di-4-ANEPPDHQ ratiometric images in
   bactch mode (whole chosen folder) obtained using a Leica microscope. The generation
   of HSB images of these GP images has been also implemented.
*/

print("\\Clear");

requires("1.44d");
closeAllImages();

// Select images folder
dir = getDirectory("Choose a Directory ");

// Inicialize choice variables
CHANNEL1 = newArray("none","ch00","ch01","ch02","ch03","ch04","ch05");
CHANNEL2 = newArray("ch00","ch01","ch02","ch03","ch04","ch05");
QUESTION = newArray("Yes","No");
Lut_Dir = getDirectory("luts");
lut = getFileList(Lut_Dir);

// Choose image channels and threshold value
Dialog.create("GP analysis parameters");
Dialog.addChoice("Acquisition ordered channel:  ", CHANNEL2, "ch00");
Dialog.addChoice("Acquisition disordered channel:  ", CHANNEL2, "ch01");
Dialog.addNumber("Lower Threshold Value for GP the mask:  ", 15);
Dialog.addChoice("Scale color for GP images:  ", lut, "Rainbow RGB.lut");
Dialog.addMessage("\n");
Dialog.addChoice("Immunofluorescence channel:  ", CHANNEL1, "none");
Dialog.addNumber("Lower Threshold Value for the IF mask:  ", 50);
Dialog.addMessage("\n");
Dialog.addNumber("G factor (1 if unknown):  ", 1);
Dialog.addChoice("Do you want to use G factor for GP image calculation?",QUESTION, "No");
Dialog.addMessage("\n");
Dialog.addChoice("Do you want to generate HSB images?",QUESTION, "Yes");
Dialog.addMessage("\n");
Dialog.show();

// Feeding variables from dialog choices
chA = Dialog.getChoice();
chB = Dialog.getChoice();
t = Dialog.getNumber();
lut = Dialog.getChoice();
chC = Dialog.getChoice();
tc = Dialog.getNumber();
Gf = Dialog.getNumber();
Ques = Dialog.getChoice();
HSB = Dialog.getChoice();

lut1= substring(lut,0,lengthOf(lut)-4);

time0 = getTime();
setBatchMode(true);

// Folder management
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
if (hour<10) {hours = "0"+hour;}
else {hours=hour;}
if (minute<10) {minutes = "0"+minute;}
else {minutes=minute;}
if (month<10) {months = "0"+(month+1);}
else {months=month+1;}
if (dayOfMonth<10) {dayOfMonths = "0"+dayOfMonth;}
else {dayOfMonths=dayOfMonth;}

var Folder = 0;
results_Dir = dir + "Results "+year+"-"+months+"-"+dayOfMonths+" "+hours+"h"+minutes+ File.separator;
File.makeDirectory(results_Dir);

ordered_images_Dir = results_Dir + "Ordered Images" + File.separator;
Folder = ordered_images_Dir;
newFolder();
disordered_images_Dir = results_Dir + "Disordered Images" + File.separator;
Folder = disordered_images_Dir;
newFolder();
GP_images_Dir = results_Dir + "GP images" + File.separator;
Folder = GP_images_Dir;
newFolder();
histogramGP_Dir = GP_images_Dir + "Histograms" + File.separator;
Folder = histogramGP_Dir;
newFolder();
rawGP_images_Dir = results_Dir + "raw GP images" + File.separator;
Folder = rawGP_images_Dir;
newFolder();

IF_images_Dir = results_Dir + "Immunofluorescence Images" + File.separator;
GP_IF_images_Dir = results_Dir + "GP-IF images" + File.separator;
histogramIF_Dir = GP_IF_images_Dir + "Histograms" + File.separator;
HSB_Dir = results_Dir + "HSB images" + File.separator;

// Open and save images in 32bit format 
listDir = getFileList(dir);
var s = 0;
for (i = 0; i < listDir.length; i++) {
	if (endsWith(listDir[i], chA+".tif")) {
		open(dir + listDir[i]);
		prepareImage();
		saveAs("tiff", ordered_images_Dir + substring(s,0,lengthOf(s)-9) + "_chA_32bits");
		close(); }
	if (endsWith(listDir[i], chB+".tif")) {
		open(dir + listDir[i]);
		prepareImage();
		Gf1=Gf;
		if (Ques=="Yes") {
			run("Multiply...","value=" + Gf);
			Gf1=1; }
		saveAs("tiff", disordered_images_Dir + substring(s,0,lengthOf(s)-9) + "_chB_32bits");
		close(); }
}

// GP and GPc Arrays 
GP=newArray(256);
for (i=0; i<256; i++) {
	GP[i]=((i-127)/127); }

GPc=newArray(256);
for (i=0; i<256; i++) {
	GPc[i]=-(1+GP[i]+Gf1*GP[i]-Gf1)/(-1-GP[i]+Gf1*GP[i]-Gf1); }
		
// GP analysis
listOrd = getFileList(ordered_images_Dir);
listDisord = getFileList(disordered_images_Dir);

var histoDir=0;
for (h = 0, j = h; h < listDisord.length; h++, j++) {
	Name=newArray(listOrd.length);
	open(ordered_images_Dir+listOrd[h]);
	name = getTitle;
	Name[j] = substring(name,0,lengthOf(name)-15);
	rename("Image_1a.tif");
	run("Duplicate...","title=Image_1b.tif");
	open(disordered_images_Dir+listDisord[j]);
	rename("Image_2a.tif");
	run("Duplicate...","title=Image_2b.tif");
	imageCalculator("Substract create 32-bit", "Image_1a.tif", "Image_2a.tif");
	rename("Image_Subs.tif");
	imageCalculator("Add create 32-bit", "Image_1b.tif", "Image_2b.tif");
	rename("Image_Add.tif");
	imageCalculator("Divide create 32-bit", "Image_Subs.tif", "Image_Add.tif");
	saveAs("tiff", rawGP_images_Dir + Name[j] + "_preGP");
	rename("Image_preGP.tif");
	setMinAndMax(-1.0000, 1.0000);
	call("ij.ImagePlus.setDefault16bitRange", 0);
	selectImage("Image_Add.tif");
	setThreshold(t*2,510);
	run("Convert to Mask");
	run("Subtract...", "value=254");
	rename("Image_1bit.tif");
	imageCalculator("Multiply create", "Image_1bit.tif", "Image_preGP.tif");
	run(lut1);
	saveAs("tiff", GP_images_Dir + Name[j] + "_GP");

// Histogram generation
	histoDir=histogramGP_Dir;
	HistogramGeneration();
}

// If there is not fluorescent immunostaining
if (chC == "none") {
// HSB image generation
	if (HSB=="Yes") {
		HSBgeneration();
	}
// Print information
	listGP = getFileList(GP_images_Dir);
	time1 = getTime();
	TOTALtime = (time1-time0)/1000;
	printInfo();
}

// If there is fluorescent immunostaining
else {
// Folder management
	temp_Dir = getDirectory("temp");

	temp_IFmask_Dir = temp_Dir + "IF mask" + File.separator;
	Folder = temp_IFmask_Dir;
	newFolder();

	File.makeDirectory(IF_images_Dir);
	File.makeDirectory(GP_IF_images_Dir);
	File.makeDirectory(histogramIF_Dir);

// Create a 1-bit mask from Immunofluorescence image and save for each image
	listDir = getFileList(dir);
	Name=newArray(listDir.length);
	for (i = 0; i < listDir.length; i++) {
		if (endsWith(listDir[i], chC + ".tif")) {
			open(dir + listDir[i]);
			name = getTitle;
			Name[i]=substring(name,0,lengthOf(name)-9);
			run("32-bit");
			saveAs("tiff",IF_images_Dir + Name[i] + "_IF_32bits");
			run("8-bit");
			setThreshold(tc,255);
			run("Convert to Mask");
			run("Divide...","value=255");
			saveAs("tiff", temp_IFmask_Dir + Name[i] + "_GP-IF_1bit");
			close();
		}
	}

// GP-IF image
	listTempIF = getFileList(temp_IFmask_Dir);
//	listL = getFileList(GP_images_Dir);
	Name=newArray(listTempIF.length);
	for (h = 0, j = h; h < listTempIF.length; h++, j++) {
		open(temp_IFmask_Dir+listTempIF[h]);
		k = getTitle;
		Name[j]=substring(k,0,lengthOf(k)-15);
		open(GP_images_Dir + Name[j] + "_GP.tif");
		l = getTitle;
		run("8-bit");
		imageCalculator("Multiply create", k, l);
		run(lut1);
		saveAs("tiff", GP_IF_images_Dir + Name[j] + "_GP-IF");

// Histogram and Normal Distribution
		histoDir=histogramIF_Dir;
		HistogramGeneration();
	}

// Folder management
	listDel = getFileList(temp_IFmask_Dir);
	for (i = 0; i < listDel.length; i++) {
		File.delete(temp_IFmask_Dir+listDel[i]); }
	File.delete(temp_IFmask_Dir);

// HSB image generation
	if (HSB=="Yes") {
		HSBgeneration();
	}

// Print information
	listGP = getFileList(GP_images_Dir);
	time2 = getTime();
	TOTALtime = (time2-time0)/1000;
	printInfo();
}


///////////////// FUNCTIONS ////////////////////

function closeAllImages() {				// This function closes all images
	while (nImages>0) {
		selectImage(nImages);
		close(); }
}

function newFolder() {					// This function creates a folder, removing any existing file in a folder with the same name
	File.makeDirectory(Folder);
	listFolder = getFileList(Folder);
	for (i = 0; i < listFolder.length; i++) {
		File.delete(Folder+listFolder[i]); }
}

function prepareImage () {				// This funcion prepares each image for the analysis
	s=getTitle;
	run("8-bit");
	run("Grays");
	run("32-bit");
	return s;
}

function HistogramGeneration () {		// This funcion obtains the intensity frequency histogram of the images, makes it smoother,
										// calculates the normal average distribution and also includes the GP value (and GP value
										// corrected by the Gfactor) corresponding to each intensity. An MS Excel file is generated
										// with all this data
	Int=newArray(256);
	Cou=newArray(256);
	Smo=newArray(256);
	NAvDist=newArray(256);
	nBins = 256;

	getHistogram(values, counts, nBins);
	close();
	while (nImages>0) {
		selectImage(nImages);
		close(); }

	for (u=0; u<nBins; u++) {
		Int[u]=u;
		Cou[u]=counts[u];
		if (u<=1) {
			Smo[u]=0; }
		else if (u==255) {
			Smo[u]=0; }
		else {
			Smo[u]=(counts[u-1]+counts[u]+counts[u+1])/3;}
	}
	Array.getStatistics(Cou,min,max,mean,stdDev);
	Sa=(mean*256)-counts[0]-counts[255];
	d=File.open(histoDir + Name[j]+"_Histogram.xls");
	print(d, "Intensity	Counts	Smooth	Norm Av Dist	GP	GP corrected");
	for (k=0; k<256; k++) {
		NAvDist[k]=100*Smo[k]/Sa;
		print(d, Int[k]+"	"+Cou[k]+"	"+Smo[k]+"	"+NAvDist[k]+"	"+GP[k]+"	"+GPc[k]);
	}
	File.close(d);
}

function HSBgeneration() {				// This function generates the HSB image of the GP images as explained in the paper
	closeAllImages();
	setBatchMode(false);

// Select images folder
	listRAW = getFileList(rawGP_images_Dir);
	
	BRIGHTNESS = newArray("Order channel","Disorder channel","IF channel");
	BRIGHTNESSnoIF = newArray("Order channel","Disorder channel");

	HSB_Dir = results_Dir + "HSB images" + File.separator;
	Folder = HSB_Dir;
	newFolder();

	Lut_Dir = getDirectory("luts");
	lut = getFileList(Lut_Dir);

	Dialog.create("Choose images and LUT");
	Dialog.addChoice("Select Hue (GP) image: ", listRAW, "none");
	if (chC == "none") {
		Dialog.addChoice("Brightness folder: ", BRIGHTNESSnoIF, "Order channel");}
	else {
		Dialog.addChoice("Brightness folder: ", BRIGHTNESS, "Order channel");}
	Dialog.addMessage("\n");
	Dialog.addChoice("Select color LUT: ",lut, "Rainbow RGB.lut");
	Dialog.addMessage("\n");
	Dialog.addChoice("Do you want to convert the whole folder into HSB images?", QUESTION, "Yes");
	Dialog.show();

// Feeding variables from dialog choices
	H = Dialog.getChoice();
	BRIGHT = Dialog.getChoice();
	Lut = Dialog.getChoice();
	WholeDir=Dialog.getChoice();

	if (BRIGHT=="Order channel") {
		mark = "_chA_32bits.tif";
		brightness_Dir = results_Dir + "Ordered Images" + File.separator;
	}
	else if (BRIGHT=="Disorder channel") {
		mark = "_chB_32bits.tif";
		brightness_Dir = results_Dir + "Disordered Images" + File.separator;
	}
	else {
		mark = "_IF_32bits.tif";
		brightness_Dir = results_Dir + "Immunofluorescence Images" + File.separator;
	}

	run("Set Measurements...", "  min limit display redirect=None decimal=5");
	index2=indexOf(Lut,".lut");
	L=substring(Lut,0,index2);

	open(rawGP_images_Dir+H);
	name=getTitle();
	Name=substring(name,0,lengthOf(name)-10);
	run(L);
	rename("Hue");

	run("Brightness/Contrast...");
	waitForUser("set min & max","set min & max");
	getMinAndMax(min,max);
	time0 = getTime();

	open(brightness_Dir+Name+mark);
	Bri=getTitle();
	run("Enhance Contrast", "saturated=0.5 normalize");

	selectWindow("Hue");
	run("RGB Color");
	run("Split Channels");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (red)");
	rename("bR");
	run("8-bit");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (green)");
	rename("bG");
	run("8-bit");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (blue)");
	rename("bB");
	run("8-bit");

	run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
	saveAs("tiff", HSB_Dir + Name + "_HSB");

	closeAllImages();

// HSB whole folder processing
	if (WholeDir == "Yes") {
		for (j = 0; j < listRAW.length; j++) {
			open(rawGP_images_Dir+listRAW[j]);
			name1=getTitle();
			Name1=substring(name1,0,lengthOf(name1)-10);
			rename("Hue");
			run(L);
			setMinAndMax(min,max);

			open(brightness_Dir+Name1+mark);
			run("Enhance Contrast", "saturated=0.5 normalize");
			Bri=getTitle;

			selectWindow("Hue");
			run("RGB Color");
			run("Split Channels");

			imageCalculator("Multiply create 32-bit", Bri,"Hue (red)");
			rename("bR");
			run("8-bit");

			imageCalculator("Multiply create 32-bit", Bri,"Hue (green)");
			rename("bG");
			run("8-bit");

			imageCalculator("Multiply create 32-bit", Bri,"Hue (blue)");
			rename("bB");
			run("8-bit");

			run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
			saveAs("tiff", HSB_Dir + Name1 + "_HSB");

			closeAllImages();
		}
	}

// Make HSB LUT bar
	newImage("Hue", "8-bit Ramp", 256, 256, 1);
	run(L);
	run("Duplicate...", "title=Brightness");
	run("Rotate 90 Degrees Left");
	run("32-bit");
	selectWindow("Hue");
	run("RGB Color");
	run("Split Channels");
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (red)");
	rename("bR");
	run("8-bit");
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (green)");
	rename("bG");
	run("8-bit");
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (blue)");
	rename("bB");
	run("8-bit");
	run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
	rename("HSB LUT");
	selectWindow("Hue (red)");
	close();
	selectWindow("Hue (green)");
	close();
	selectWindow("Hue (blue)");
	close();
	selectWindow("Brightness");
	close();

	selectWindow("HSB LUT");
	run("Rotate 90 Degrees Left");
	run("Size...", "width=32 height=256 interpolation=None");


// Copy LUT bar to the image
	run("Copy");
	newImage("Panel LUT", "RGB White", 70, 264, 1);
	setBatchMode(false);
	makeRectangle(4,4,32,256);
	run("Paste");
	run("Colors...", "foreground=black background=black selection=yellow");
	run("Line Width...", "line="+2);
	run("Draw");
	run("Select None");
	setFont("Arial", 12);
	setColor(0, 0, 0);
	drawString(d2s(max,2),39,15);
	drawString(d2s(min,2),39,264);
	run("Select None");

	if (WholeDir == "Yes") {
		saveAs("tiff", HSB_Dir + "Lut bar");
		close();
	}
	else {
		saveAs("tiff", HSB_Dir + Name + "_Lut bar");
		open(HSB_Dir + Name + "_HSB.tif");
	}
	closeAllImages();
}

function printInfo () {					// This function prints and saves a summary of the results of the macro
	print("\\Clear");
	print("GP images analysis macro");
	print("   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms");
	print("   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus");
	print("   Nature Protocols 2011");
	print("   version July 2011");
	print("\n");
	print("Ordered channel: " + chA);
	print("Disordered channel: " + chB);
	print("   Lower threshold value: " + t);
	if (chC != "none") {
		print("Channel IF: " + chC);
		print("   Lower threshold value IF mask: " + tc); }
	print("");
	print("GP images (" + (listGP.length-1) + ") saved at:");
	print("  " + GP_images_Dir);
	print("G factor: " + Gf);
	print("Scale color for GP images: " + lut1);
	print("");
	if (chC != "none") {
		listGP2 = getFileList(GP_IF_images_Dir);
		print("GP-IF images (" + (listGP2.length-1) + ") saved at:");
		print("  " + GP_IF_images_Dir);
		print(""); }
	if (HSB=="Yes") {
		print("HSB images saved at:");
		print("  " + HSB_Dir);
		print("");
	}
	print("");
	print("Execution time: " + d2s(TOTALtime,2) + " seg.");
	print("");
	print("Analysis date: "+DayNames[dayOfWeek]+", "+dayOfMonth+" "+MonthNames[month]+" "+year+" - "+hours+":"+minutes);
	print("ImageJ version " + getVersion());

	selectWindow("Log");
	saveAs("Text", results_Dir + "Log GP images generation "+year+"-"+months+"-"+dayOfMonths+" "+hours+"h"+minutes);

	setBatchMode("exit and display");	
}
