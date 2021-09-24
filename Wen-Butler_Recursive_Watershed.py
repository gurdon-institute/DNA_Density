# Recursive watershed segmentation within a tolerance range to segment diffuse nuclei during nuclear membrane breakdown
# 	by Richard Butler, Gurdon Institute Imaging Facility

import math as maths

from ij import IJ, ImagePlus
from ij.gui import Roi, TextRoi, ShapeRoi, Overlay
from ij.process import ImageStatistics, Blitter, ImageProcessor, AutoThresholder, FloodFiller
from ij.plugin.filter import ThresholdToSelection, MaximumFinder, EDM
from ij.measure import ResultsTable



def fillHoles(mask):
	width = mask.getWidth()
	height = mask.getHeight()
	ff = FloodFiller(mask)
	mask.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if mask.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if mask.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if mask.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if mask.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	n = width*height
	for i in range(n):
		if mask.get(i)==127:
		    mask.set(i, 0)
		else:
		    mask.set(i, 255)


def watershed(ip, tol):
	floatEdm = EDM().makeFloatEDM(ip, 0, False)
	maxIp = MaximumFinder().findMaxima(floatEdm, tol, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
	if (maxIp != None):
		ip.copyBits(maxIp, 0, 0, Blitter.AND)


def watershedROI(roi, tol):
	bounds = roi.getBounds()
	mask = roi.getMask()
	watershed(mask, tol)
	parts = getRois(mask)
	for part in parts:
		offset = part.getBounds()
		part.setLocation(bounds.x+offset.x, bounds.y+offset.y)
	return parts


def getRois(mask):
	rois = []
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = ThresholdToSelection().convert(mask)
	rois = ShapeRoi(composite).getRois()
	return rois


# segment using recursive watershed with decreasing tolerance until area is acceptable
def segment(roi, tol):
	segmented = []
	area = roi.getStatistics().area * cal.pixelWidth * cal.pixelHeight
	if area < minA:
		return []
	elif area < segmentA:
		return [roi]
	else:
		subRois = watershedROI(roi, tol)
		for sub in subRois:
			subarea = sub.getStatistics().area * cal.pixelWidth * cal.pixelHeight
			if tol >= minTol and subarea > segmentA:
				segmented.extend( segment(sub, tol-1) )
			elif subarea >= minA:
				segmented.append(sub)
	return segmented


imp = IJ.getImage()
title = imp.getTitle()
cal = imp.getCalibration()
ol = Overlay()

sigma = 0.4/cal.pixelWidth
minA = 15.0		#µm²
segmentA = 50.0
minTol = 4
startTol = 12

c1ip = imp.getStack().getProcessor(1)
c2ip = imp.getStack().getProcessor(2)

ip = c1ip.duplicate()
ip.blurGaussian(sigma)
hist = ip.getHistogram(256)
stats = ip.getStatistics()
thresh = AutoThresholder().getThreshold( AutoThresholder.Method.Huang, hist )
thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
ip.threshold(int(thresh))
mask = ip.convertToByte(False)
fillHoles(mask)

nuclei = []
for roi in getRois(mask):
	nuclei.extend( segment(roi, startTol) )

rt = ResultsTable.getResultsTable()
for i,nuc in enumerate(nuclei):
	ol.add(nuc)
	row = rt.getCounter()
	rt.setValue("Image", row, imp.getTitle())
	rt.setValue("Nucleus", row, i)
	area = nuc.getStatistics().area
	rt.setValue("Area ("+u"\u00b5"+"m"+u"\u00b2", row, area*cal.pixelWidth*cal.pixelHeight)
	perim = nuc.getLength()
	circ = 4*maths.pi*(area/(perim*perim))
	rt.setValue("Circularity", row, circ)
	c1ip.setRoi(nuc)
	c1Stats = c1ip.getStatistics()
	c2ip.setRoi(nuc)
	c2Stats = c2ip.getStatistics()
	rt.setValue("C1 Mean", row, c1Stats.mean)
	rt.setValue("C1 StdDev", row, c1Stats.stdDev)
	rt.setValue("C2 Mean", row, c2Stats.mean)
	rt.setValue("C2 StdDev", row, c2Stats.stdDev)

imp.setOverlay(ol)
rt.show("Results")
