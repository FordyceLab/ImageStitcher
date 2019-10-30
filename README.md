# ImageStitcher
!["One Piece at a Time"](/resources/one_piece_at_a_time.png)

## Overview
ImageStitcher is a set of classes and functions to stitch (concatenate tiled arrays), flat-field correct, and background-subtract image rasters. It is primarily intended for use with images generated via Fordyce lab RunPack, but is also capable of stitching arbitrary multi-dimensional Micromanager-generated .ome.tif stacks (TODO). Basic logging is provided for convenience.

### Architecture
- **StitchingSettings**: High-level class for management of general stiching settings and flat-field correction images/parameters
- **RasterParams**: A struct-like class for raster acquisition parameter handling
- **Raster**: An image raster class. Tile images stored in memory or as references.
- **RasterSet**: An un-ordered collection class for arbitrary Raster objects.
- **RasterGroup**: A general raster superclass to handle associated (grouped) rasters.
- **KineticImaging**: A RasterGroup subclass specific to kinetic imaging (time indexed)
- **FileHandler**: A generic file handling class for basic folder traversal and (future) handling of various image storage schema (for MicroManager support, for example)
- **BackgroundImages**: A simple class for handling "Background" images and execute background subtraction

- **Utils (To Add)**: Basic Jupyter Notebooks for FF correction parameter calculation and EXIF metadata extraction

## Tasks
- [ ] Clean up and add utils 
