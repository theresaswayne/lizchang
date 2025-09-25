#@ File	(label = "Input directory", style = "directory") srcFile
#@ File	(label = "Output directory", style = "directory") dstFile
#@ String  (label = "File extension", value=".zip") ext
#@ String  (label = "File name contains", value = "RoiSet") containString
#@ Integer  (label = "ROIs to skip at beginning", value = "1") skipRois
#@ Integer (label = "# Pixels to dilate nucleus before creating cytoplasm", value = "3") dilate
#@ boolean (label = "Keep directory structure when saving", value = true) keepDirectories

# Cytoplasm_ROI.py
# Given an ROIset of a prescribed format, containing nucleus and whole cell, 
# adds new ROIs for cytoplasm only
# measures area of cell and cytoplasm in pixels
# ------- EXPECTED ROISet FORMAT ----------
# optional: 1st ROI = background measurement
# 2nd and 3rd: nucleus and whole cell, respectively, for cell # 1
# 4th and 5th, etc.: nucleus and cell for subsequent cells

# Theresa Swayne, 2025
#  -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

from ij import IJ, ImagePlus, ImageStack
from ij.plugin import ZProjector
from ij.plugin.filter import RankFilters
# from ij.plugin.filter import BackgroundSubtracter
import net.imagej.ops
from net.imglib2.view import Views
from net.imglib2.img.display.imagej import ImageJFunctions as IL
from ij.process import ImageStatistics as IS
from net.imglib2.algorithm.dog import DogDetection
from ij.gui import PointRoi
from jarray import zeros
from ij.measure import ResultsTable
from math import sqrt
from java.awt import Color
from ij.plugin.frame import RoiManager
from ij.gui import GenericDialog
import os
from loci.plugins import BF
from jarray import array
from ij.measure import ResultsTable

def process(imp, srcDir, dstDir, currentDir, fileName, keepDirectories, skip, dilate, table):
	
 	#IJ.run("Close All", "")
	
	imp = IJ.getImage() # necessary to avoid RoiMgr errors


	rm = RoiManager.getInstance()
	if not rm:
	  rm = RoiManager()
	rm.reset()
	
 	IJ.log("Processing ROI set:" + fileName)
 	
 	rm.runCommand("Open", os.path.join(currentDir, fileName))
 	
 	# loop through ROIs, skipping as needed
 	
 	numRois = rm.getCount()
 	IJ.log("This set has " + str(numRois) + " ROIs")
 	
 	# check for a valid number of ROIs
 	if (numRois - skipRois) % 2 != 0:
 		IJ.log("Skipped ROI set " + fileName + " because it has an invalid number of ROIs.")
 		return

	startIndex = skipRois # indices start at 0, so if we skip one, we start at 1
 	endIndex = numRois - startIndex # end is not included in range 
 	cellCount = 1
 	for RoiIndex in range(startIndex, endIndex, 2): # if we skip the first ROI, indices will be 1, 3, etc
		IJ.log("Processing ROI index " + str(RoiIndex))
		rm.select(RoiIndex)
		IJ.run("Enlarge...", "enlarge="+str(dilate))
		#RoiEnlarger.enlarge(imp, 20)
		rm.runCommand(imp,"Update")
 		rm.rename(RoiIndex, "Nucl_" + str(cellCount))
		rm.rename(RoiIndex+1, "Cell_" + str(cellCount))
		cellRois = array([RoiIndex, RoiIndex+1], 'i')
 		rm.setSelectedIndexes(cellRois)
 		IJ.log("Selected ROIs " + str(cellRois[0]) + "," + str(cellRois[1]))
 		rm.runCommand("XOR") # requires an image to be open
 		rm.runCommand("Add") # should be added at the end of the list
 		rm.deselect() # make sure nothing else selected
 		newTotal = rm.getCount()
 		IJ.log("There are now " + str(newTotal) + " ROIs")
 		rm.rename(newTotal-1, "Cyto_" + str(cellCount))
 		cellCount = cellCount + 1
 	
 	
 	# measure area
 	numRois = rm.getCount()
 	
 	for RoiIndex in range(0, numRois):
	 		
	 	# add a line to the results table
		table.incrementCounter()
		table.addValue("Filename", fileName)
		
		rm.select(RoiIndex);
		roi = imp.getRoi()
		roiName = roi.getName()
		table.addValue("ROI name", roiName)
		rm.runCommand(imp,"Measure");
		stats = imp.getStatistics(IS.AREA)
		#IJ.log("area: %s" %(stats.area))
	
		# Add to results table
		table.addValue("ROI area",stats.area)
		rm.deselect() # make sure nothing else selected

	# save the updated ROIs and area table
	saveDir = currentDir.replace(srcDir, dstDir) if keepDirectories else dstDir
	if not os.path.exists(saveDir):
		os.makedirs(saveDir)
	IJ.log("Saving ROIs to" + saveDir)
	baseName = os.path.splitext(fileName)[0]
	rm.deselect() # make sure nothing else selected
	rm.save(os.path.join(saveDir, baseName + "_Cyto_Rois.zip"))
	# table.save(os.path.join(saveDir, baseName + "_RoiAreas.csv"))


def run():
	srcDir = srcFile.getAbsolutePath()
	dstDir = dstFile.getAbsolutePath()
	
	IJ.log("\\Clear")
	IJ.log("Processing ROIsets")
	
	imp = IJ.createImage("Dummy", "8-bit black", 2048, 2048, 1) # necessary to avoid RoiMgr errors
	imp.show()
	
	table = ResultsTable()
	
	for root, directories, filenames in os.walk(srcDir):
		filenames.sort();
	for filename in filenames:
		# Check for file extension
		if not filename.endswith(ext):
			continue
		# Check for file name pattern
		if containString not in filename:
			continue
		process(imp, srcDir, dstDir, root, filename, keepDirectories, skipRois, dilate, table)

	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.reset()
	
	table.save(os.path.join(dstDir, "RoiAreas.csv"))
	
	IJ.run("Close All", "")
	IJ.run("Clear Results")
	IJ.log("Done!")

run()

#RoiManager.getName(index)