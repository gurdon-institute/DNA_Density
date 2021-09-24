import math as maths
import re

from ij import IJ, ImagePlus
from ij.gui import Roi, TextRoi, ShapeRoi, Overlay
from ij.process import ImageStatistics, Blitter, ImageProcessor, AutoThresholder, FloodFiller
from ij.plugin.filter import ThresholdToSelection, MaximumFinder, EDM
from ij.measure import ResultsTable

from java.awt import Color, Font, Canvas


def getRoi(initialRoi, channel, method):
	ip = imp.getStack().getProcessor(channel).duplicate()
	ip.blurGaussian(sigma)
	ip.setRoi(initialRoi)
	hist = ip.getHistogram(256)
	stats = ip.getStatistics()
	thresh = AutoThresholder().getThreshold( method, hist )
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
	ip.setThreshold(thresh, 2E16, ImageProcessor.NO_LUT_UPDATE)
	roi = ThresholdToSelection().convert(ip)
	return roi

def getDensity(roi, channel):
	ip = imp.getStack().getProcessor(channel).duplicate()
	ip.findEdges()
	ip.setRoi(roi)
	stats = ip.getStatistics()
	density = (stats.mean+stats.stdDev) / maths.sqrt(stats.area)
	return density
	

imp = IJ.getImage()
title = imp.getTitle()
cal = imp.getCalibration()
sigma = 0.25/cal.pixelWidth
ol = Overlay()

cellRoi = getRoi(None, 1, AutoThresholder.Method.Minimum)	# find total signal area
#cellRoi.setStrokeColor(Color.CYAN)
#ol.add(cellRoi)
cellA = cellRoi.getStatistics().area

nLevels = 3
C1rois = [None for t in range(nLevels)]
colours = [Color(0.5,0.5,1.0), Color.GREEN, Color.YELLOW]
for t in range(nLevels):
	prev = C1rois[t-1] if t>0 else None
	roi = getRoi(prev, 1, AutoThresholder.Method.Otsu)
	if t == 0:
		roi.setStrokeColor(colours[t])
		ol.add(roi)
	C1rois[t] = roi

hcRoi = getRoi(C1rois[0], 2, AutoThresholder.Method.MaxEntropy)	#map heterochromatin in C2 inside lowest density C1 chromatin ROI
hcRoi.setStrokeColor(Color.MAGENTA)
ol.add(hcRoi)

imp.setOverlay(ol)

C1andC2rois = [ ShapeRoi(roi).and(hcRoi) for roi in C1rois]

areaC1 = [roi.getStatistics().area for roi in C1rois]
areaC2 = hcRoi.getStatistics().area
areaC1C2 = [roi.getStatistics().area for roi in C1andC2rois] 

time = re.findall("[0-9]+h_", title)[0]
if time:
	time = time[:-2]
else:
	time = "?"

rt = ResultsTable.getResultsTable()
i = 0
row = rt.getCounter()
rt.setValue("Image", row, title)
rt.setValue("Time (h)", row, time)
rt.setValue("Cell Area ("+cal.getUnit()+u"\u00b2"+")", row, cellA*cal.pixelWidth*cal.pixelHeight)
rt.setValue("Density Level", row, i)
rt.setValue("DNA Density", row, getDensity(C1rois[i], 1) )
rt.setValue("DNA Area ("+cal.getUnit()+u"\u00b2"+")", row, areaC1[i]*cal.pixelWidth*cal.pixelHeight )
rt.setValue("Proportion of Cell Area", row, areaC1[i]/cellA )
rt.setValue("Heterochromatin Area ("+cal.getUnit()+u"\u00b2"+")", row, areaC2 )
rt.setValue("Heterochromatin Proportion", row, areaC1C2[i]/areaC1[i] )
rt.setValue("DNA Density in Heterochromatin", row, getDensity(C1andC2rois[i], 1) )
	
rt.show("Results")
