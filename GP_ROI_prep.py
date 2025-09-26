#@ File	(label = "Input directory", style = "directory") srcFile
#@ File	(label = "Output directory", style = "directory") dstFile
#@ String  (label = "File extension", value=".zip") ext
#@ String  (label = "File name contains", value = "RoiSet") containString
#@ Integer  (label = "ROIs to skip at beginning", value = "1") skipRois
#@ boolean (label = "Keep directory structure when saving", value = true) keepDirectories

# GP_ROI_prep.py
# Given an ROIset of a prescribed format, containing cytoplasm and membrane ROIs, f
# renames ROIs according to cell number and ROI type
# ------- EXPECTED ROISet FORMAT ----------
# optional: 1st ROI = background measurement
# 2nd and 3rd: cytoplasm and membrane, respectively, for cell # 1
# 4th and 5th, etc.: cytoplasm and membrane for subsequent cells
# -------LIMITATIONS---------------
# for speed, generates a dummy image of 2048 x 2048, and cannot handle ROIs in larger images

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
from ij.measure import ResultsTable
from math import sqrt
from java.awt import Color
from ij.plugin.frame import RoiManager
from ij.gui import GenericDialog
import os
from loci.plugins import BF
from jarray import array, zeros


def process(imp, srcDir, dstDir, currentDir, fileName, keepDirectories, skip):
	
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
 		rm.rename(RoiIndex, "Cyto_" + str(cellCount))
		rm.rename(RoiIndex+1, "Memb_" + str(cellCount))
 		rm.deselect() # make sure nothing else selected
 		cellCount = cellCount + 1

	# remove stack position info so we can see all ROIs when Show All is selected
	rm.deselect() # make sure nothing selected so we operate on all ROIs
	rm.runCommand(imp,"Remove Channel Info");
	rm.runCommand(imp,"Remove Slice Info");
	rm.runCommand(imp,"Remove Frame Info");


	# save the updated ROIs
	saveDir = currentDir.replace(srcDir, dstDir) if keepDirectories else dstDir
	if not os.path.exists(saveDir):
		os.makedirs(saveDir)
	IJ.log("Saving ROIs to" + saveDir)
	baseName = os.path.splitext(fileName)[0]
	rm.deselect() # make sure nothing is selected
	rm.save(os.path.join(saveDir, baseName + "_Rois.zip"))


def run():
	srcDir = srcFile.getAbsolutePath()
	dstDir = dstFile.getAbsolutePath()
	
	IJ.log("\\Clear")
	IJ.log("Processing ROIsets")
	
	imp = IJ.createImage("Dummy", "8-bit black", 2048, 2048, 1) # dummy image to allow ROI operations
	imp.show()
	
	for root, directories, filenames in os.walk(srcDir):
		filenames.sort();
	for filename in filenames:
		# Check for file extension
		if not filename.endswith(ext):
			continue
		# Check for file name pattern
		if containString not in filename:
			continue
		process(imp, srcDir, dstDir, root, filename, keepDirectories, skipRois)

	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.reset()
	
	IJ.run("Close All", "")
	IJ.run("Clear Results")
	IJ.log("Done!")

run()

#RoiManager.getName(index)