# ImageStitcher [Rework, Minimal Testing]
!["One Piece at a Time"](/resources/one_piece_at_a_time.png)

## Overview
ImageStitcher is a set of classes and functions to stitch (concatenate tiled arrays), flat-field correct, and background-subtract image rasters. It is primarily intended for use with images generated via Fordyce lab RunPack, but is also capable of stitching arbitrary multi-dimensional Micromanager-generated .ome.tif stacks. Basic logging is provided for convenience.

### Architecture
- **StitchingSettings**: High-level class for management of general stiching settings and flat-field correction images/parameters
- **RasterParams**: 
- **Raster**: 
- **RasterSet**:
- **RasterGroup**: A general raster superclass to handle associated rastered images
- **KineticImaging**: A RasterGroup subclass specific to kinetic imaging
- **FileHandler**:
- **BackgroundImages**:

- **Utils (To Add)**: Basic Jupyter Notebooks for FF correction parameter calculation and EXIF metadata extraction

### Setup
#### Configuration
1. config.json: 

#### Notebook
Coming Soon!
## Tasks
- [ ] Clean up and add utils 
