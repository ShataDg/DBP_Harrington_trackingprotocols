# To run this script you need pyimagej and openjdk=8. The best way I found to run was a new enviroment
#conda https://pyimagej.readthedocs.io/en/latest/Install.html#installing-via-conda-mamba
#conda activate pyimagej

# To test your installation: python -c 'import imagej; ij = imagej.init("2.5.0"); print(ij.getVersion())'

#more information can be found at pyimagej Documentation (https://pyimagej.readthedocs.io/en/latest/Install.html#installing-via-conda-mamba)

import os
import sys
from tokenize import Double3
import imagej
import scipy.ndimage
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scyjava as sj
from jpype import JInt
from imagej._java import jc
sj.config.add_option('-Xmx6g')

import code # for debugging -- remove later

run_mode = 'headless' # modes should be "headless" or "interactive"

# initialize fiji endpoint to get plugins
ij = imagej.init('sc.fiji:fiji',mode=run_mode)
print(f"ImageJ version: {ij.getVersion()}")

# Setup the params as you seem fit
params = {
    'time': "30",  # Time interval between frames in minutes
    'outradius': 4,  # Remove outliers radius
    'outthreshold': 10,  # Remove outliers threshold
    'outwhich': "Bright",  # Remove outliers which (Bright/Dark)
    'diameter': 20.0,  # Expected object diameter
    'threshold': 0.5,  # Quality threshold
    'medianfilter': True,  # Pre-process with median filter
    'subpixel': True,  # Sub-pixel localization
    'channel': 1,  # Channel
    'quality': 1.2,  # Spot detection quality filter
    'link_distance': 15.0,  # Linking max distance; Default value is 15
    'gap_distance': 15.0,  # # Gap closing will not occur for two spots if their distance exceeds this value; Default value is 15
    'splitting_max_distance': 15.0, # Track splitting will not occur for two spots if their distance exceeds this value; Default value is 15
    'merging_max_distance': 15.0, # Track merging will not occur for two spots if their distance exceeds this value; Default value is 15
    'cutoff_percentile': 0.9, # Default value is 0.9; takes double
    'alternative_linking_cost_factor': 1.05, # Default value is 1.05; takes double
    'default_blocking_value': float('inf'),
    'gap_frame': JInt(20),  # Gap-close max frame gap
    'input_dir': '/Volumes/BAGGINS/TRD2_Harrington_Data/DBPHarrington/Results/input_snc_complete',  # Output folder
    'output_dir': '/Volumes/BAGGINS/TRD2_Harrington_Data/DBPHarrington/Results/output_snc_complete',  # Input folder (manually set)
}
print(params['input_dir'])

# legacy ImageJ resources
DirectoryChooser = sj.jimport('ij.io.DirectoryChooser')
GenericDialog = sj.jimport('ij.gui.GenericDialog')
ResultsTable = sj.jimport('ij.measure.ResultsTable')
TextWindow = sj.jimport('ij.text.TextWindow')

# ImageJ2/Fiij resources
Model = sj.jimport('fiji.plugin.trackmate.Model')
Settings = sj.jimport('fiji.plugin.trackmate.Settings')
TrackMate = sj.jimport('fiji.plugin.trackmate.TrackMate')
DogDetectorFactory = sj.jimport('fiji.plugin.trackmate.detection.DogDetectorFactory')
SparseLAPTrackerFactory = sj.jimport('fiji.plugin.trackmate.tracking.jaqaman.SparseLAPTrackerFactory')
#LAPUtils = sj.jimport('fiji.plugin.trackmate.tracking.LAPUtils')
LabelImgExporter = sj.jimport('fiji.plugin.trackmate.action.LabelImgExporter')
label_id_painting = sj.jimport('fiji.plugin.trackmate.action.LabelImgExporter.LabelIdPainting')
DisplaySettingsIO = sj.jimport('fiji.plugin.trackmate.gui.displaysettings.DisplaySettingsIO')
FeatureFilter = sj.jimport('fiji.plugin.trackmate.features.FeatureFilter')
Table = sj.jimport('org.scijava.table.Table')

# Java resrouces
Double = sj.jimport('java.lang.Double')

# Indicating which params to use depending upon the mode (headless or interactive) selected
def setup():
    setup_params = {}
    if run_mode !="headless": 
        dc = DirectoryChooser("Choose a folder")
        path = dc.getDirectory()
        if path is None:
            print("User Canceled")
            return
        else:
            gd = GenericDialog(path)
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
            gd.addStringField("Output folder", "/home/user/output/", 60)
            gd.showDialog()
        if gd.wasCanceled():
            print("Setup canceled")
            return
        else:
            setup_params['time'] = str(gd.getNextNumber())
            setup_params['outradius'] = str(gd.getNextNumber())
            setup_params['outthreshold'] = str(gd.getNextNumber())
            setup_params['outwhich'] = str(gd.getNextChoice())
            setup_params['diameter'] = float(gd.getNextNumber())
            setup_params['threshold'] = float(gd.getNextNumber())
            setup_params['medianfilter'] = gd.getNextBoolean() 
            setup_params['subpixel'] = gd.getNextBoolean()
            setup_params['channel'] = JInt(gd.getNextNumber())
            setup_params['quality'] = float(gd.getNextNumber())
            setup_params['link_distance'] = float(gd.getNextNumber())
            setup_params['gap_distance'] = float(gd.getNextNumber())
            setup_params['gap_frame'] = JInt(gd.getNextNumber())
            setup_params['output_dir'] = str(gd.getNextString())
            setup_params['input_dir'] = str(path)
            return setup_params
    else:
        #code to set the params
        return params


def detection(params) -> dict:
    detection_results = {}
    file_type = '.tif'
    files = os.listdir(params['input_dir'])
    for f in files:
        if f[-4:] == file_type:
            name = f.split(".")[0]
            imp = ij.IJ.openImage(os.path.join(params['input_dir'], f))
            imp.show()
            ij.IJ.selectWindow(name + file_type)
            
            # get number of frames/slices
            # Assuming 3D data, the first two dimensions in ImageJ are always
            # X and Y, thus the last dimension will be frames/slices
            for i in range(len(imp.dims)):
                if (imp.dims[i] != 'X') and (imp.dims[i] != 'Y'):
                    frames = imp.shape[i]
            ij.IJ.run(imp, "Properties...", "channels=1 slices=1 frames=" + str(frames) + " pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000 frame=[" + params['time'] +" min]")
            ij.IJ.run(imp, "Despeckle", "stack")
            ij.IJ.run("Despeckle", "stack")
            ij.IJ.run("Despeckle", "stack")
            ij.IJ.run(imp, "Remove Outliers...", "radius=" + str(params['outradius']) + " threshold=" + str(params['outthreshold']) + " which=" + str(params['outwhich']) +" stack")
            detection_results['raw_image'] = ij.WindowManager.getImage(name + '.tif').duplicate() # save raw/processed image

            model = Model()
            settings = Settings(imp)
            radius = params['diameter'] / 2
            settings.detectorFactory = DogDetectorFactory()
            settings.detectorSettings = {
                'DO_SUBPIXEL_LOCALIZATION': params['subpixel'],
                'RADIUS': radius,
                'TARGET_CHANNEL': params['channel'],
                'THRESHOLD': params['threshold'],
                'DO_MEDIAN_FILTERING': params['medianfilter']
            }
            filter = FeatureFilter('QUALITY', params['quality'], True)
            settings.addSpotFilter(filter)
            settings.trackerFactory = SparseLAPTrackerFactory()
            #settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
            settings.trackerSettings['LINKING_MAX_DISTANCE'] = params['link_distance']
            settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = params['gap_distance']
            settings.trackerSettings['ALLOW_GAP_CLOSING'] = True
            settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
            settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = params['splitting_max_distance']
            settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
            settings.trackerSettings['MERGING_MAX_DISTANCE'] = params['merging_max_distance']
            settings.trackerSettings['CUTOFF_PERCENTILE'] = params['cutoff_percentile']
            settings.trackerSettings['ALTERNATIVE_LINKING_COST_FACTOR'] = params['alternative_linking_cost_factor']
            settings.trackerSettings['BLOCKING_VALUE'] = params['default_blocking_value']
            settings.trackerSettings['MAX_FRAME_GAP'] = params['gap_frame']
            tm = TrackMate(model, settings)

            # check TrackMate input and process
            if not tm.checkInput():
                sys.exit(str(tm.getErrorMessage()))

            tm.process()
            #use_spots_as_id = True
            export_spots = True
            export_tracks = True
            label_id_painting = LabelImgExporter.LabelIdPainting.LABEL_IS_INDEX
            label_img = LabelImgExporter.createLabelImagePlus(tm, export_spots, export_tracks, label_id_painting)
            detection_results['label_image'] = label_img.duplicate() # save label image
            detection_results['TRACK_ID'] = None
            label_img.show()
            ij.IJ.saveAs("Tiff", os.path.join(params['output_dir'], name + '_lblImg.tif'))
            ij.IJ.run("Close")

            export_spots = False
            export_tracks = True
            label_id_painting = LabelImgExporter.LabelIdPainting.LABEL_IS_INDEX
            spot_img = LabelImgExporter.createLabelImagePlus(tm, export_spots, export_tracks, label_id_painting)
            spot_img.show()
            ij.IJ.saveAs("Tiff", os.path.join(params['output_dir'], name + '_SpotImg.tif'))
            ij.IJ.run(imp, "Merge Channels...", "c4=" + name + ".tif c6=" + name + "_SpotImg.tif keep ignore")
            ij.IJ.run(imp, "RGB Color", "frames")
            ij.IJ.saveAs("Tiff",os.path.join(params['output_dir'], name + '_merge.tif'))
            ij.IJ.run("Close")

            fm = model.getFeatureModel()
            rt = ij.WindowManager.getWindow("TrackMate Results")

            # table stuff
            if (rt == None) or not (isinstance(rt, TextWindow)):
                table = ResultsTable()
            else:
                table = rt.getTextPanel().getOrCreateResultsTable()
            table.reset

            for id in model.getTrackModel().trackIDs(True):
                v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
                model.getLogger().log('')
                model.getLogger().log('Track ' + str(id) + ': mean velocity = ' + str(v) + ' ' + str(model.getSpaceUnits()) + '/' + str(model.getTimeUnits()))
                track = model.getTrackModel().trackSpots(id)
                for spot in track:
                    sid = spot.ID()
                    x = spot.getFeature('POSITION_X')
                    y = spot.getFeature('POSITION_Y')
                    t = spot.getFeature('FRAME')
                    q = spot.getFeature('QUALITY')
                    model.getLogger().log('\tspot ID = ' + str(sid) + ': x=' + str(x) + ', y=' + str(y) + ', t=' + str(t) + ', q=' + str(q))
                    table.incrementCounter()
                    table.addValue('TRACK_ID', Double(id))
                    table.addValue('SPOT_ID', sid)
                    table.addValue('POSITION_X', x)
                    table.addValue('POSITION_Y', y)
                    table.addValue('FRAME', t)
                    table.addValue('QUALITY', q)

            if run_mode != "headless":
                # below code is applicable only for interactive
                table.show("TrackMate Results")
                ij.WindowManager.getWindow("TrackMate Results").rename(name + '.csv')
                ij.IJ.saveAs("Measurements", os.path.join(params['output_dir'], name + '.csv'))
                ij.IJ.selectWindow(name + '.csv')
                ij.IJ.run("Close")

            else:
                table.show("TrackMate Results")
                #ij.IJ.saveAs("Measurements", os.path.join(params['output_dir'], name + '.csv'))
                table.save(os.path.join(params['output_dir'], name + '.csv'))

            detection_results['raw_image'].show()
            ij.IJ.selectWindow('DUP_' + name + '.tif')
            ij.IJ.saveAs("Tiff",os.path.join(params['output_dir'], name + '_processed.tif'))
            detection_results['tracks'] = ij.py.from_java(ij.convert().convert(table, Table)) # convert ResultsTable to pandas dataframe

            return detection_results

def process_tracks(
    df: pd.DataFrame,
    label_image,
    raw_image,
    max_event_frame: int = 500,
    min_event_frames_for_good: int = 60,
    event_ext: int = 0,
    pixel_size: int = 11
    ):
    # convert ImagePlus to numpy.ndarray
    if isinstance(label_image, jc.ImagePlus):
        label_image = ij.py.from_java(label_image).data
    if isinstance(raw_image, jc.ImagePlus):
        raw_image = ij.py.from_java(raw_image).data
    
    length = raw_image.shape[0]
    c = int(raw_image.shape[1] / 2)
    raw_image_left = raw_image[:, :, :c]
    raw_image_right = raw_image[:, :, c:]

    # sort the input dataframe for event ROI based on event number and frame number
    # check if TRACK_ID is none - raise ValueType Error
    df = df.where(df.notnull(), None)
    df = df.sort_values(by=['TRACK_ID', 'FRAME'])
    df['TRACK_ID'] = df['TRACK_ID'] +  1
    df['FRAME'] = df['FRAME'] + 1

    # filter events around cell movements
    tracks  = df.groupby('FRAME').agg(**{'TRACK_COUNT': pd.NamedAgg(column='TRACK_ID', aggfunc='nunique'), 'TRACK_LIST': pd.NamedAgg(column='TRACK_ID', aggfunc='unique')})
    multi_tracks = tracks.drop(tracks[tracks.TRACK_COUNT <= max_event_frame].index)
    multi_tracks = multi_tracks.explode('TRACK_LIST')
    multi_tracks_list = multi_tracks.TRACK_LIST.unique()
    bad_tracks = multi_tracks_list.tolist()

    # save the x, y position at the start and end of the track for each event
    track_dict = {}
    for index, row in df.iterrows():
        if row['TRACK_ID'] != None:
            if row['TRACK_ID'] not in bad_tracks:
                track_id = int(row['TRACK_ID'])
                if track_id not in list(track_dict.keys()):
                    track_dict[track_id] = {}
                    track_dict[track_id]['start_frame'] = track_dict[track_id]['stop_frame'] = row['FRAME']
                    track_dict[track_id]['start_x'] = track_dict[track_id]['stop_x'] = row['POSITION_X']
                    track_dict[track_id]['start_y'] = track_dict[track_id]['stop_y'] = row['POSITION_Y']
                else:
                    track_dict[track_id]['stop_frame'] = row['FRAME']
                    track_dict[track_id]['stop_x'] = row['POSITION_X']
                    track_dict[track_id]['stop_y'] = row['POSITION_Y']
    track_list = list(track_dict.keys())
    
    # filter tracks that are just 1 frame
    for t in track_list:
        if int(track_dict[t]['stop_frame']) - int(track_dict[t]['start_frame']) < 1:
            del track_dict[t]
    track_list = list(track_dict.keys())

    # add ROI in the same x, y position at the start and end of the track event_ext number of frames
    for t in track_list:
        start = int(track_dict[t]['start_frame']) - event_ext
        stop = int(track_dict[t]['stop_frame']) + event_ext + 1
        if start <= 0:
            start = 1
        else:
            start = start
        if stop > length:
            stop = length
        else:
            stop = stop

        for frame in range(start, stop, 1):
            if t not in label_image[frame]:
                label_image[frame] - np.where(label_image[frame - 1] == t, t, label_image[frame])

    # dilate ROI to pixel_size of event
    for frame in range(length):
        a = label_image[frame]
        b = scipy.ndimage.maximum_filter(a, pixel_size)
        b[a != 0] = a[a != 0]
        label_image[frame] = b

    # save measurement of event pre filtering
    pre_filter_df = pd.DataFrame(index=track_list)
    for frame in range(length):
        pre_filter_df[frame] = scipy.ndimage.mean(raw_image[frame], labels=label_image[frame], index=track_list)
    pre_filter_df = pre_filter_df.T
    pre_filter_df['timepoint'] = pre_filter_df.index
    pre_filter_df.dropna(axis='columns', how='all', inplace=True)
    pre_filter_df.to_csv('pre_filter.csv')

    # filter good cells/tracks -- here you can choose the min number of frames for good tracks
    good_tracks = []
    for t in track_list:
        start = int(track_dict[t]['start_frame'])
        stop = int(track_dict[t]['stop_frame'])
        l = stop - start
        if l > min_event_frames_for_good:
            good_tracks.append(t)
    
    post_filter_df = pd.DataFrame(index=track_list)
    post_filter_df = post_filter_df[post_filter_df.index.isin(good_tracks)]

    # measure the mean intensity of good tracks
    for frame in range(length):
        post_filter_df[frame] = scipy.ndimage.mean(raw_image[frame], labels=label_image[frame], index=good_tracks)
    post_filter_df = post_filter_df.T
    post_filter_df.dropna(axis='columns', how='all', inplace=True)
    post_filter_df.to_csv('post_filter.csv')

    # measure the mean intensity of the entire, left and right side of the image
    y0 = []
    y1 = []
    y2 = []
    for i in range(raw_image.shape[0]):
        y0.append(raw_image_left[i, :, :].mean())
        y1.append(raw_image_right[i, :, :].mean())
        y2.append(raw_image[i, :, :].mean())

    y0_df = pd.DataFrame(y0)
    y1_df = pd.DataFrame(y1)
    y2_df = pd.DataFrame(y2)

    post_filter_df['left'] = y0_df
    post_filter_df['right'] = y1_df
    post_filter_df['whole'] = y2_df

    # save mean intensity of good tracks
    #post_filter_df.to_csv(os.path.join(output_folder, 'complete_' + movie + '.csv'))
    post_filter_df.plot.line(legend=True)
    plt.show()

    return

if __name__ == "__main__":
    params = setup()
    print(params)
    detection_results = detection(params)
    process_tracks(detection_results['tracks'], detection_results['label_image'], detection_results['raw_image'])

    # return REPL
    code.interact(local=locals())
    
