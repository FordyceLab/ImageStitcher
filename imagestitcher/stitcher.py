# title             : stitcher.py
# description       : General-Purpose Rastered Image Stitching
# authors           : Daniel Mokhtari
# credits           : Craig Markin, Polly Fordyce
# date              : 20180520
# version update    : 20180520
# version           : 0.1.0
# usage             : With permission from DM
# python_version    : 3.6


# General Python
import os
import re
import pathlib
import warnings
import logging
from copy import deepcopy
from datetime import datetime

# Scientific Data Structures and Plotting
import numpy as np
import pandas as pd
from matplotlib import pyplot as pl
from skimage import io, transform, external


channels = {'1pbp', '2bf', '3dapi', '4egfp', '5cy5'}

class StitchingSettings:
    ffPaths = None
    ffParams = None
    ffImages = None
    rasterPattern =  (True, False) # Raster origin. topleft = (0, 1), topright = (1, 0), 
    tileDim = None

    def __init__(self, ffPaths = None, ffParams = None, tileDim = 1024, setupNum = 1):
        """
        StitchingSettings for general stitching parameters

        Arguments:
            (dict) ffPaths:
            (ffParams) ffParams:
            (int) tileDim:
            (int) setupNum:

        Returns:
            None
        """
        StitchingSettings.ffPaths = self.defineFFPaths(ffPaths)
        StitchingSettings.ffImages = self.ffImages = self.readFFImages()
        StitchingSettings.ffParams = self.ffParams = ffParams
        StitchingSettings.tileDim = self.tileDimensions = tileDim
        if setupNum != 2: 
            StitchingSettings.rasterPattern = (True, False)
        else:
            StitchingSettings.rasterPattern = (False, True)
        self.initializeLogger()


    def readFFImages(self):
        """
        Loads FF images

        Arguments:
            None

        Returns:
            None
        """

        result = {}
        for channel, path in self.ffPaths.items():
            result[channel] = io.imread(path)
        return result


    def defineFFPaths(self, ffChannelsPaths):
        """
        Error checking for FF paths dictionary

        Arguments:
            (dict) ffChannelsPaths: dictionary of channels:ffpath

        Returns:
            (dict) Dictionary of channel:ffpath

        """
        allPaths = {}
        for channel, path in ffChannelsPaths.items():
            if os.path.isfile(path): 
                allPaths[channel] = path
            else: 
                note = 'Error: The Flat-Field Image for {} at {} does not exist'
                raise ValueError(note.format(channel, path))
        return allPaths


    def showFFImages(self, vmin = 0, vmax = 65535):
        """
        Displays loaded FF images
        
        Arguments:
            (int) vmin: intensity minimum
            (int) vmax: intensity maximum

        Returns;
            None

        """
        for channel, image in StitchingSettings.ffImages.items():
            fig = pl.imshow(image, cmap = 'gray', vmin = vmin, vmax = vmax)
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('RFU')
            pl.suptitle('{} FF Correction Image'.format(channel))
            pl.show()


    def initializeLogger(self):
        """
        Initializes the logger

        Arguments:
            None

        Return:
            None
        """
        logfmt = '%(asctime)s %(levelname)-8s %(message)s'
        logging.basicConfig(format=logfmt, level=logging.DEBUG, datefmt='%y-%m-%d %H:%M:%S')


class RasterParams:
    def __init__(self, root, source, overlap, exposure, channel, setup, 
        groupFeature = 0, rotation = 0, autoFF = True):
        """
        Parameters describing a single image raster.

        Arguments:
            (str) root: root path
            (str) source: image source ('ipnb' | 'mm')
            (float) overlap: overlap fraction (e.g., 0.1)
            (int) exposure: exposure time (ms)
            (str) channel: imaging channel ('1pbp' | '2bf' | '3dapi' | '4egfp' | '5cy5')
            (int) setup: setup number (1 | 2 | 3 | 4)
            (int | float | string) groupFeature: feature value for RasterGroup
            (float) rotation: pre-stitch rotation to perform (%)
            (bool) autoFF: flag to execute FF correction on stitch, if possible

        Returns:
            None

        """
        self.root = root
        self.parent = pathlib.Path(root).parents[0]
        self.source = source #ipnb or mm
        self.size = StitchingSettings.tileDim
        self.dims = None
        self.overlap = overlap
        self.exposure = exposure
        self.channel = channel
        self.setup = setup
        self.rotation = rotation
        self.acquiOri = StitchingSettings.rasterPattern
        self.groupFeature = groupFeature
        self.autoFF = autoFF



class FileHandler:
    def __init__(self):
        return


    def parseKineticFolder(self, paramTemplate, path):
        """
        I-Python Specific: Parses a folder containing time-series imaging

        Arguments:
            (RasterParams) paramTemplate: a template RasterParams object
            (str) path: root path containg kinetic imaging folders
        
        Returns:
            (KineticImaging) A KineticImaging object of timecourse imaging

        """
        root = pathlib.Path(path)
        subfolders = [i for i in root.iterdir() if not 'StitchedImages' in i.name]
        def getTime(pathName):
            chipAcqTime = re.match('^\d{8}_\d{6}|^\d{8}-\d{6}', pathName)
            acquiTime =  datetime.strptime(chipAcqTime[0], '%Y%m%d-%H%M%S')
            return acquiTime
        
        getRasterFolders = lambda f: (f, getTime(f.name), self.parseFolder(f))
        rasterfolders = [getRasterFolders(f) for f in subfolders]

        rasters = set
        for subfolder, time, imageSet in rasterfolders:
            pt = deepcopy(paramTemplate)
            pt.groupFeature = time
            pt.root = subfolder
            dims, exposuresImages = self.readImages(imageSet, 'paths') #Grab the raster dimensions
            rasters = rasters.union(self.makeRasters(pt, exposuresImages, dims).rasters)
        return KineticImaging(path, rasters)


    def parseSingleFolder(self, paramTemplate, path):
        """
        I-Python Specific: Parses a folder containing a single channel of rastered images
        (can contain multiple exposure times) and generates a set of rasters (a RasterSet)
        ready to stitch and export.

        Arguments:
            (RasterParams) paramTemplate: a template RasterParams object
            (str) path: root path containg images
        
        Returns:
            (RasterSet) A RasterSet containing a bag of Rasters

        """

        imageSet = self.parseFolder(path)
        dims, exposuresImages = self.readImages(imageSet, 'paths')
        return self.makeRasters(paramTemplate, exposuresImages, dims)


    def parseFolder(self, path):
        """
        I-Python Specific: Parses a folder containing a single channel of rastered images
        (can contain multiple exposure times) and generates a DataFrame of indexed 
        tile paths. For using in Raster generation.

        Arguments:
            (str) path: root path containg images
        
        Returns:
            (pd.DataFrame) A DataFrame of indexed tiled image paths

        """
        root = pathlib.Path(path)
        paths = [pathlib.Path(i) for i in root.iterdir()]
        fileNames = [path.parts[-1] for path in paths]

        match = '\d{3}_\d{3}_\d{1,4}'
        indexFiles = lambda f: tuple(re.search(match, f)[0].split('_'))

        raster = {'paths': paths, 'files': fileNames}
        image = pd.DataFrame(raster)

        image['indexes'] = image.files.apply(indexFiles)
        image['x'] = image.indexes.apply(lambda i: int(i[0]))
        image['y'] = image.indexes.apply(lambda i: int(i[1]))
        image['exp'] = image.indexes.apply(lambda i: int(i[2]))
        return image.sort_values(['x', 'y']).reset_index(drop = True).set_index('exp')


    def makeRasters(self, paramTemplate, exposuresImages, dims):
        """
        I-Python Specific: Generates a RasterSet from a DataFrame of indexed rastered image paths.
        
        Arguments:
            (RasterParams) paramTemplate: a template RasterParams object
            (dict) exposuresImages: a dictionary of exposure:list of rastered images
            (tuple) dims: raster tiling dimensions (not image dimensions)
        
        Returns:
            (RasterSet) A collection of Rasters 

        """
        rasters = []
        for exp, images in exposuresImages.items():
            pc = deepcopy(paramTemplate)
            pc.exposure = exp
            pc.dims = dims
            rasters.append(FlatRaster(images, pc))
        return RasterSet(set(rasters))
    

    def readImages(self, df, feature):
        """
        I-Python Specific: Extracts a raster worth of image path, in proper raster order, from a folder 
        record generated by parseFolder
        
        Arguments:
            (pd.DataFrame) df: DataFrame containing indexed tile paths
            (str) feature: name of DataFrame column containg path (typ. 'paths')

        Returns:
            (tuple) A tuple of raster tiling dimensions and ordered rastered image paths


        """
        readImage = lambda i: io.imread(i)

        # Each folder can have multiple exposure times, but one channel
        exposures = set(df.index)
        imageRefs = {}
        for exp in exposures:
            workingRaster = df.loc[exp]
            xdims = workingRaster['x'].drop_duplicates().size
            ydims = workingRaster['x'].drop_duplicates().size
            dims = (xdims, ydims)
            imageRefs[exp] = [i for i in workingRaster[feature].tolist()]
        return (dims, imageRefs) #exposure time: list of image paths




class Raster:
    def __init__(self, imageRefs, params):
        """
        A collection of images (or references) for one rastered image at one set of acquisition
        parameters

        Arguments:
            (list) imageRefs: an ordered list of rastered image paths
            (RasterParams) params: rastered image parameters

        Returns:
            None


        """
        self.imageRefs = imageRefs
        self.params = params

    def applyFF(self):
        """
        Applies flat-field correction to fetched images

        Arguments:
            None

        Returns:
            None
        
        """
        channel = self.params.channel
        exposure = self.params.exposure
        ffImage = StitchingSettings.ffImages[channel]
        ffParams = StitchingSettings.ffParams[channel][exposure]
        ffbval = ffParams[0]
        ffscale = ffParams[1]
        def ffSubtract(i, ffi, ffbval, ffscale):
            ffresult = np.subtract(i, ffbval)/np.subtract(ffi, ffbval)*ffscale
            result = np.clip(ffresult, 0, 65535).astype('uint16')
            return result
        return [ffSubtract(i, ffImage, ffbval, ffscale) for i in self.images]


    def stitch(self, method = 'cut'):
        """
        Wrapper for image stitching method selection.
        TODO: Implement 'overlap' method

        Arguments:
            (str) method: stitch method ('cut' | 'overlap')
        
        Returns:
            (np.ndarray) A stitched image array

        """

        self.fetchImages()
        if method == 'cut':
            return self.cutStitch()
        elif method == 'overlap':
            return self.overlapStitch()
        else:
            raise ValueError('Invalid stitch method. Valid methods are "cut" and "overlap"')


    def cutStitch(self):
        """
        Stitches a raster via the 'cut' method. Trims borders according to image overlap and 
        concatenates along both axes. If RasterParameters.autoFF is True, performs 
        flat-field correction prior to stitching.

        Arguments:
            None

        Returns:
            (np.ndarray) A stitched image array

        """

        imsize = self.params.size
        margin = int(imsize*self.params.overlap/2) #Edge dim to trim
        retained = imsize-2*margin #Remaining tile dim
        border = slice(margin,-margin)

        tiles = self.images
        if (self.params.autoFF 
            and self.params.channel in StitchingSettings.ffParams 
            and self.params.exposure in StitchingSettings.ffParams[self.params.channel]):
            tiles = self.ffCorrectedImages = self.applyFF()

        trimmedTiles = [tile[border, border] for tile in tiles] #Trim
        tileArray = np.asarray(trimmedTiles) #Make ndarray
        arrangedTiles = np.reshape(tileArray, (self.params.dims[0], self.params.dims[0], retained, retained))
        if self.params.acquiOri[0]: #If origin on right, flip horizontally (view returned)
            arrangedTiles = np.flip(arrangedTiles, 0)
        if self.params.acquiOri[1]: #If origin on bottom, flip vertically (view returned)
            arrangedTiles = np.flip(arrangedTiles, 1)
        rowsStitched = np.concatenate(arrangedTiles, axis = 2) #Stitch rows
        fullStitched = np.concatenate(rowsStitched, axis = 0) #Stitch cols
        return fullStitched


    def overlapStitch(self):
        """
        #TODO: re-implement overlap stitching method

        """
        
        return


    def exportStitch(self, method = 'cut', outPathName = 'StitchedImages'):
        """
        Perform stitching and export raster.

        Arguments:
            (str) method: stitch method ('cut' | 'overlap')
            (str) outPathName: Name of folder to house stitched raster. Typically 'StitchedImages'

        Returns:
            None

        """
        stitchedRaster = self.stitch(method = method)
        
        features = [self.params.exposure, self.params.channel, self.params.groupFeature]
        rasterName = 'StitchedImage_{}_{}_{}.tif'.format(*features)
        
        stitchDir = pathlib.Path(os.path.join(self.params.parent, outPathName))
        stitchDir.mkdir(exist_ok = True)
        outDir = os.path.join(stitchDir, rasterName)
        external.tifffile.imsave(outDir, stitchedRaster)
        logging.info('Stitching Complete')




class FlatRaster(Raster):
    def fetchImages(self):
        """
        Fetches (loads into memory) and rotates images (if indicated by raster parameters)

        Arguments:
            None

        Returns:
            None

        """
        r = self.params.rotation
    
        rotationParams = {'resize': False, 'clip': True, 'preserve_range': True}
        rotateImage = lambda i: transform.rotate(i, r, **rotationParams).astype('uint16')
        
        readImage = lambda i: io.imread(i)
        images = [readImage(i) for i in self.imageRefs]
        
        if r:
            images = [rotateImage(i) for i in images]
        
        self.images = images



class StackedRaster(Raster):
    def __init__(self, imageRefs, stackIndices, params):
        """

        Arguments:
        
        Returns:
            None

        """
        super().__init__(imageRefs, stackIndices)
        self.stackIndices = stackIndices


    def fetchImages(self):
        """
        Fetches (loads into memory) and rotates images (if indicated by raster parameters)

        Arguments:
            None

        Returns:
            None

        """
        r = self.params.rotation
    
        rotationParams = {'resize': False, 'clip': True, 'preserve_range': True}
        rotateImage = lambda i: transform.rotate(i, r, **rotationParams).astype('uint16')
        
        readImage = lambda i: io.MultiImage(i[1])i[0] #Second elem is ref, first is stack no.
        images = [readImage(i) for i in zip(self.stackIndices self.imageRefs)]
        
        if r:
            images = [rotateImage(i) for i in images]
        
        self.images = images






class RasterSet:
    def __init__(self, rasters):
        """
        A set of rasters

        """
        self.rasters = rasters

    def exportStitchAll(self, **args):
        """
        Basic stitching and export of raster set
        """
        while self.rasters:
            self.rasters.pop().exportStitch(**args)



class RasterGroup:
    def __init__(self, root, rasters):
        """


        """
        # Dict of feature (time) 
        self.root = root
        self.rasters = rasters #list of rasters


    def add(self):
        #add raster to group
        return


    def stitch(self):
        while self.raster:
            self.raster.pop().stitch()



class KineticImaging(RasterGroup):
    def order(self):
        """
        Orders the set of rasters as a dictionary of time:raster entries

        Arguments:
            None

        ReturnsL:
            None

        """

        sortedTimes = sorted([raster.params.groupFeature for raster in self.rasters])
        self.referencedTimes = [t - min(sortedTimes) for t in sortedTimes]
        self.orderedRasters = {raster.params.groupFeature-min(sortedTimes):raster for raster in self.rasters}


    def exportStitch(self, method = 'cut', outPathName = 'StitchedImages2'):
        """
        Stitches and exports each of the rasters, appending the raster time onto the filename

        Arguments:
            (str) method: stitch method ('cut' | 'overlap')
            (outPathname) outPathName: Name of folder to house stitched raster. Typically 'StitchedImages'

        Returns:
            None

        """
        for dt, raster in self.orderedRasters.items():
            time = dt.total_seconds()
            stitchedRaster = raster.stitch(method = method)
            
            features = [raster.params.exposure, raster.params.channel, int(time)]
            rasterName = 'StitchedImage_{}_{}_{}.tif'.format(*features)
            

            stitchDir = pathlib.Path(os.path.join(self.root, outPathName))
            stitchDir.mkdir(exist_ok = True)
            outDir = os.path.join(stitchDir, rasterName)
            external.tifffile.imsave(outDir, stitchedRaster)

            mp = {'t': time, 'ch':raster.params.channel, 'ex':raster.params.exposure}
            logging.info('Stitch Saved: (Time: {t}, Ch: {ch}, Ex: {ex})'.format(**mp))
        logging.info('Kinetic Stitches Complete')



class BackgroundImages:
    def __init__(self):
        """
        A simple class for background subtraction.
        Stores a collection of background images and permits scripted background subtraction

        Arguments:
            None

        Returns:
            None

        """

        self.images = {}


    def add(self, path, index, channel, exposure):
        """
        Adds a background image, described by an index described by (index, channel, exposure)
        where index is typically device name (i.e., 'd1')

        Arguments:
            (str) path: Background image path
            (str) index: arbitrary index (typically device name (i.e., 'd1'))
            (str) channel: Background image channel
            (int) exposure: Background image exposure time (ms)

        Returns:
            None

        """

        self.images[(index, channel, exposure)] = io.imread(path)


    def remove(self, index, channel, exposure):
        """
        Removes a background image.
        Arguments:
            (str) index: arbitrary index (typically device name (i.e., 'd1'))
            (str) channel: Background image channel
            (int) exposure: Background image exposure time (ms)

        Returns:
            None

        """

        del(self.images[(index, channel, exposure)])


    def substractBackground(self, targetImagePath, targetindex, targetchannel, targetexposure, prefix = 'BGSubstracted_'):
        """
        Subtracts a background image from a target image of matching index, channel, and exposure.
        Resulting images are prefixed with an optional prefix, typically 'BGSubtracted_'

        Arguments:
            (str) targetImagePath: path of image to background subtract
            (str) targetindex: arbitrary index to match to backgrounds (typ. device name (i.e., 'd1'))
            (str) targetchannel: target image channel
            (int) targetexposure:: target image exposure
            (str) prefix: resulting image filename prefix

        Returns:
            None


        """
        logging.info('Background Substrating Stitch: (Ch: {ch}, Ex: {ex})'.format(**mp))
        imgDir = pathlib.Path(targetImagePath)
        target = io.imread(imgDir)

        bgImage = self.image[(targetindex, targetchannel, targetexposure)]
        bgsub = np.subtract(target, bgImage)
        bgsubclipped = np.clip(bgsub, 0, 65535).astype('uint16') 
        
        outPath = os.path.join(tp.parents[0], '{}{}'.format(prefix, imgDir.name))
        external.tifffile.imsave(outPath, bgsubclipped)

        mp = {'ch': targetchannel, 'ex': targetexposure}
        logging.info('Background Subtraction Complete')


    def walkAndBGSubtract(self, path, index, channel):
        """
        Walks a directory structure, find images to background subtract, and executes subtraction

        Arguments:
            (str) path: path from hwere to walk
            (str) index: arbitrary index to select background image
            (str) channel: channel to select background image
        
        Returns:
            None

        """

        correctable = lambda f: channel in f and 'StitchedImage' in f and not 'BGSubtracted' in f
        parse = lambda f: tuple(f.split('.')[0].split('_')[-3:])
        
        for root, dirs, files in os.walk(path):
            if 'StitchedImages' in root:
                toCorrect = {parse(f):f for f in files if correctable(f)}
                for params, file in toCorrect.items():
                    exposure, channel, feature = params
                    self.substractBackground(file, index, channel, exposure)


########## Standard scripts ##########
def stitchImage(path, params):
    """
    Stitch and export a single raster at the path with given parameters (exposure is overwritten)

    """

    mp = {'ch': params.channel, 'ex': params.exposure, 'o':params.overlap, 'r':params.rotation}
    startmessage = 'Stitching Folder (Ch: {ch}, Exp: {ex}, Overlap: {o}, Rot: {r})'.format(**mp)
    logging.info(startmessage)
    fh = FileHandler()
    raster = fh.parseSingleFolder(params, path)
    raster.exportStitchAll()


def stitchKinetics(path, params):
    """
    Stitch and export a timecourse of rasters at the path with given parameters

    """

    startmessage = 'Starting Kinetic Stitch'
    logging.info(startmessage)
    fh = FileHandler()
    k = fh.parseKineticFolder(params, path) #Returns KineticImaging
    k.order()
    k.exportStitch()


def walkAndStitch(path, params, stitchtype = 'kinetic'):
    """
    Walk a directory structure of identical types (all single acquisitions or all kinetic)
    and stitch all images beneath.

    Arguments:
        (RasterParams) params: raster parameters. Exposure time ignored.
        (str) stitchtype: type of rasters contained ('single' | 'kinetic')

    Returns:
        None

    """

    global channels
    for root, dirs, files in os.walk(path):
        basename = os.path.basename(root)
        tidydirs = [direct for direct in dirs if 'StitchedImages' not in direct]
        if basename in channels:
            channel = os.path.basename(root)
            if stitchtype == 'kinetic':
                newParams = deepcopy(params)
                newParams.channel = basename
                stitchKinetics(root, newParams)
            elif stitchtype == 'single':
                for direct in tidydirs:
                    newParams = deepcopy(params)
                    newParams.channel = basename
                    newParams.root = direct
                    stitchImage(direct, newParams)
                