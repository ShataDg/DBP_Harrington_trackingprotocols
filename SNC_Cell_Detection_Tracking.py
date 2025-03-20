# This script removes noise and perfoms cell detection and tracking in movies from SNC regions
# This code was tested in ImageJ 1.54p, Java 1.8.0_322 using TrackMate v7.14.0

# What this code does:
	# The code will pre process the movies (using frames and time to the stack the images), and will remove noise using Despeckle and Remove Outliers.
	# Next, it will detect cells using Trackmate DOG detector, track cells using Simple LAP tracker along the movie, save the movies (in .tif format) and the export the data in a .csv file.

# When you run the code:
	# Input: a window will pop-up asking you to point it to the folder with the images
# This code assumes that all the movies are in tif format and inside the same folder
# The movie in this case is a image stack containing all the images of a particular time-series.
	# Set the following parameters
		# Time interval between frames in minutes (default value = 30)
		# Remove outliers radius (default value = 4)
		# Remove outliers threshold (default value = 10)
		# Remove outliers which (default value = Bright)
		# Expected object diameter: (default value = 20) Value of the cell diameter in pixels
		# Quality Threshold: (default value = 0.5) discard spurious spots 
		# Pre-process with median filter: (default value = True) this will apply a 3x3 median filter prior to any processing. This can help dealing with images that have a marked salt and pepper noise which generates spurious spots.
		# Sub-pixel localization: (default value = True)
		# Channel (default value = 1) If you're using a multichannel image you can choose the channel of interest here.
		# Initial thresholding, based on Quality (default value = 1.2) Set a threshold value based on the quality feature to limit the number of spots before calculating other features, filtering spots can save computational time, spots discarded here will not be saved.
		# Linking max distance (default value = 15) distance in pixels, This value limits the spatial search range for candidate matching spots. If you know the maximal displacement of your object between two frame, you can put this value here
		# Gap-closing max distance (default value = 15) distance in pixels, Two track segments will not be bridged if the last spot of the first segment is further away than the first spot of the second segment. 
		# Gap-close max frame gap (default value = 20)frame interval, sets the maximal frame interval between two spots to be bridged. 
		# Output folder": provide the path to the output folder

	# Output
	# The following files should be saved to the output folder 
		# 1. moviename_lblmg.tif, (movie with the centroid of the cells detected, this file will be used to detect the cells where the measure will be made)
		# 2. moviename_merge.tif, (movie containig the processed image merged with the cells detected, just for checking)
		# 3. moviename_SpotImg.tif, (movie containg just the masks of the cells detected)
		# 4. moviename_processed.tif, (movie processed from remove noise step, this file will be used to measure the mean intensity of each cell along time)
		# 5. moviename.csv (file containing the track information, this file will be used as a track reference for the intensity measurements)


#import libraries

import sys, os
import math

from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser
from ij.gui import WaitForUserDialog, GenericDialog
from ij.plugin import ImageCalculator, PlugIn
from ij.measure import Measurements, ResultsTable
from ij.plugin.frame import RoiManager

from fiji.plugin.trackmate import Model, Settings, TrackMate, SelectionModel, Logger, Spot, SpotCollection
from fiji.plugin.trackmate.detection import DogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.action import ExportAllSpotsStatsAction
from fiji.plugin.trackmate.action import LabelImgExporter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider
from fiji.plugin.trackmate.providers import EdgeAnalyzerProvider
from fiji.plugin.trackmate.providers import TrackAnalyzerProvider
from fiji.plugin.trackmate.action.LabelImgExporter.LabelIdPainting import LABEL_IS_INDEX
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter


#File type
file_type = ".tif" #file type if using .tif (you can change the code here if you're using a different file format) TIFF, GIF, JPEG, “raw”, FITS, PGM or PNG format.

# Choose the input and output folders and Trackmate parameters
dc=DirectoryChooser("Choose a folder")
dir=dc.getDirectory() #folder were all time-lapse images are saved
if dir is None:
    print("User Canceled")

# User input parameters to run the code
else:
  gd=GenericDialog(dir)
  gd.addNumericField("Time interval between frames in minutes", 30,0)
  gd.addNumericField("Remove outliers radius", 4,0)
  gd.addNumericField("Remove outliers threshold", 10,0)
  gd.addChoice("Remove outliers which", ["Bright", "Dark"], "Bright")
  gd.addNumericField("Expected object diameter",20,2)
  gd.addNumericField("Quality Threshold",0.5,1)
  gd.addCheckbox("Pre-process with median filter:", True)
  gd.addCheckbox("Sub-pixel localization:", True)
  gd.addNumericField("Channel",1,0)
  gd.addNumericField("Spot detection Quality filter",1.2,1)
  gd.addNumericField("Linking max distance",15,1)
  gd.addNumericField("Gap-closing max distance",15,1)
  gd.addNumericField("Gap-close max frame gap",20,0)
  gd.addStringField("Output folder", "/UsersDesktop/output", 60)
  gd.showDialog()  
  if gd.wasCanceled():
    print("User canceled dialog!")
  else:
 
    
    time = str(gd.getNextNumber())
    outradius = str(gd.getNextNumber())
    outthreshold = str(gd.getNextNumber())
    outwhich = str(gd.getNextChoice())
    diameter = float(gd.getNextNumber())
    threshold = float(gd.getNextNumber()) #threshold for the event detection larger values will result in less events detected, bright events being detected,
                #smaller values will result in more event but this might include false events
    medianfilter = gd.getNextBoolean()
    subpixel = gd.getNextBoolean()
    channel = int(gd.getNextNumber()) #trackmate can work on images with more then one channel if that is the case select the channel were the events are being detected
    quality = float(gd.getNextNumber()) #trackmate can filter spots based on quality, choose a value
    link_distance = int(gd.getNextNumber()) #this set the distance that an event can whitin a frame for static events this value should be set close to 0
    gap_distance = int(gd.getNextNumber()) #this set the distance that an event can move from frame to frame, for static events this value should be set close to 0
    gap_frame = int(gd.getNextNumber()) #number of frames that can have no event detected between two event and still be considered the same event
    output_dir=gd.getNextString() #folder were files created will be saved


ic =  ImageCalculator();
reload(sys)
sys.setdefaultencoding('utf-8')

# Movie processing:
#1. read all the files in the folder
#2. assign frames and time to the stack
#3. Remove some noise from images
#4. Create a Trackmate model/settings
#5. Export movies
#6. Save the tracks

#1. read all the files in the folder
movies = os.listdir(dir)
print (movies)
for file in movies:
	if file[-4:]== file_type:
		name = file.split(".")[0]
		
		#open files
		imp = IJ.openImage(os.path.join(dir,file));
		imp.show()
		IJ.selectWindow(name+'.tif');
		
		#get dimensions and set the number of frames
		n=imp.getDimensions();
		frames=str(max(n[3],n[4]));
		
		#2. stack pre processing to add frames and time and remove some noise from images
		IJ.run(imp, "Properties...", "channels=1 slices=1 frames="+frames+" pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000 frame=["+time+" min]");
		IJ.run(imp, "Despeckle", "stack");
		IJ.run("Despeckle", "stack");
		IJ.run("Despeckle", "stack");
		IJ.run(imp, "Remove Outliers...", "radius="+outradius+" threshold="+outthreshold+" which="+outwhich+" stack");
		IJ.saveAs("Tiff",os.path.join(output_dir,name+'_processed.tif'));

		#Create a trackmate model/settings
		model = Model()
		settings = Settings(imp)
		radius = float(diameter)/2
		
		settings.detectorFactory = DogDetectorFactory()
		settings.detectorSettings = {
   				'DO_SUBPIXEL_LOCALIZATION' : subpixel,
    			'RADIUS' : float(radius),
    			'TARGET_CHANNEL' : int(channel),
    			'THRESHOLD' : float(threshold),
    			'DO_MEDIAN_FILTERING' : medianfilter,
				}
		filter1 = FeatureFilter('QUALITY', quality, True)
		settings.addSpotFilter(filter1)	
			
		settings.trackerFactory = SparseLAPTrackerFactory()
		settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
		settings.trackerSettings['LINKING_MAX_DISTANCE'] = float(link_distance)
		settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = float(gap_distance)
		settings.trackerSettings['MAX_FRAME_GAP'] = int(gap_frame)
		trackmate = TrackMate(model, settings)
		
		ok = trackmate.checkInput()
		if not ok:
			sys.exit(str(trackmate.getErrorMessage()))
		
		ok = trackmate.process()
		if not ok:
			sys.exit(str(trackmate.getErrorMessage()))
		
		selectionModel = SelectionModel(model)
		
		#Export movies
		exportSpotsAsDots = True
		exportTracksOnly = True
		lblImg = LabelImgExporter.createLabelImagePlus(trackmate, exportSpotsAsDots, exportTracksOnly, LABEL_IS_INDEX)
		lblImg.show()
		IJ.saveAs("Tiff",os.path.join(output_dir,name+'_lblImg.tif'));
		IJ.run ("Close");
		
		exportSpotsAsDots = False
		exportTracksOnly = True
		SpotImg = LabelImgExporter.createLabelImagePlus(trackmate, exportSpotsAsDots, exportTracksOnly, LABEL_IS_INDEX)
		SpotImg.show()
		IJ.saveAs("Tiff",os.path.join(output_dir,name+'_SpotImg.tif'));
		IJ.run(imp, "Merge Channels...", "c4="+name+"_processed.tif c6="+name+"_SpotImg.tif keep ignore");
		IJ.run(imp, "RGB Color", "frames");
		IJ.saveAs("Tiff",os.path.join(output_dir,name+'_merge.tif'));
		IJ.run ("Close");
		IJ.run ("Close");
		
		# Save the detected tracks for each detected ROI
		ds = DisplaySettingsIO.readUserDefault()
		
		fm = model.getFeatureModel()
		
		rt_exist = WindowManager.getWindow("TrackMate Results")
		
		if rt_exist==None or not isinstance(rt_exist, TextWindow):
			table= ResultsTable()
		else:
			table = rt_exist.getTextPanel().getOrCreateResultsTable()
		table.reset
		
		for id in model.getTrackModel().trackIDs(True):
			v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
			model.getLogger().log('')
			model.getLogger().log('Track ' + str(id) + ': mean velocity = ' + str(v) + ' ' + model.getSpaceUnits() + '/' + model.getTimeUnits())
			track = model.getTrackModel().trackSpots(id)
		
			for spot in track:
				sid = spot.ID() 
				x=spot.getFeature('POSITION_X')
				y=spot.getFeature('POSITION_Y')
				t=spot.getFeature('FRAME')
				q=spot.getFeature('QUALITY')
				snr=spot.getFeature('SNR')
				model.getLogger().log('\tspot ID = ' + str(sid) + ': x='+str(x)+', y='+str(y)+', t='+str(t)+', q='+str(q))
				table.incrementCounter()
				table.addValue("TRACK_ID", id)
				table.addValue("SPOT_ID", sid)
				table.addValue("POSITION_X", x)
				table.addValue("POSITION_Y", y)
				table.addValue("FRAME", t)
				table.addValue("QUALITY", q)
		table.show("TrackMate Results") 
		IJ.selectWindow('TrackMate Results');
		IJ.saveAs("Measurements",os.path.join(output_dir,name+'.csv'));
		IJ.run("Close All");
		IJ.run("Close All");
		
IJ.run("Close All");
print("Finished");