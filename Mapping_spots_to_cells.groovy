// VisiumHDAnalysis - Quantify VisiumHD H&E Slide 
// 
// By: Ofra Golani
//
// Workflow
// ============================================================
// This code runs on a QuPath script editor. The code produces detections (spots) and annotations (cells) with measurements associating spots to cells.
// Specifically, the code does the following steps:
// - Segment Tissue region using pixel classifier
// - Segment Cells within Tissue using StarDist with Expansion
// - Convert Cells to Annotations 
// - Load VisiumHD spots as detections
// - Update Hirarchy - to associate spots to cells 
// - (Optional) Run Pixel Classifier that marks the tissue under the VisiumHD grid (WholeTissue)


import qupath.ext.stardist.StarDist2D
import groovy.time.*
import java.io.BufferedReader;
import java.io.FileReader;
import qupath.lib.objects.PathAnnotationObject;
import qupath.lib.objects.PathDetectionObject;
import qupath.lib.objects.PathTileObject;
import qupath.lib.roi.ROIs
import qupath.lib.roi.RectangleROI

// ===================  Workflow control Parameters  ====================================

def segmentTissue             = 0 // 
def segmentCells              = 0 // 
def loadSpots                 = 0
def associateSpotsToCells     = 0
def runPixelClassifierForSpot = 0
def runPixelClassifierForCell = 0 // MAKE SURE THAT AddMeasurementsToCells IS ALSO 1
def AddMeasurementsToCells    = 0
def exportCellLabels          = 0 // keep 0
def exportSpotLabels          = 0 // keep 0
def exportCellsAsGeoJson      = 1
def saveResultTable           = 1 // export result table to Tab-separated txt file

var cellClassName = "Cell"
//var nucClassName  = "Nuc"
var spotClassName = "Spot"

def resultsSubFolder = 'results' // subfolder for table 

def wholeTissueClass = "WholeTissue" 

def WholeTissueClassifier = "WholeTissue"  //M2// name of WholeTissue pixel classifier
//ef WholeTissueClassifier = "WholeTissue_M1"  //M1 // name of WholeTissue pixel classifier
// def WholeTissueClassifier = "WholeTissue"  //M2// name of WholeTissue pixel classifier

def WholeTissue_MinSize =  25          // Minimal WholeTissue connected-component size        
def WholeTissue_MinHoleSize = 50        // minmal hole size to keep when creating WholeTissue regions, samller holes are filled 

// ===================  Cell Segmentation parameters  ====================================
var pathModel = 'A:/shared/QuPathScriptsAndProtocols/QuPath_StarDistModels/he_heavy_augment.pb'
// Stardist parameters 
var clear_existing_detections = false
var param_threshold    = 0.5 //threshold for detection. All cells segmented by StarDist will have a detection probability associated with it, where higher values indicate more certain detections. Floating point, range is 0 to 1. Default 0.5
var normalize_low_pct  = 1   //lower limit for normalization. Set to 0 to disable
var normalize_high_pct = 99  // upper limit for normalization. Set to 100 to disable.
var param_expansion    = 6 //5   //size of cell expansion in pixels. Default is 10.
var param_tilesize     = 1024 //size of tile in pixels for processing. Must be a multiple of 16. Lower values may solve any memory-related errors, but can take longer to process. Default is 1024.


// ===================  Pixel Classifier Parameters (optional)   ====================================
def PixelClassifier = "xxx"

// ===================  Visium HD Spots parameters  ====================================

//double spot_diameter_fullres = 6.241827204164965 // M2 ; // obtain this value from the file scalefactors_json.json
//double spot_diameter_fullres = 6.223977607892008 //M1
double spot_diameter_fullres = 6.211844235377622 //M6

 
def csvfile = "//outs/binned_outputs/square_002um/spatial/tissue_positions.csv" // M6 - SUBSET IMAGE

// ======================================================================================================
// ===================  Code Begins - Dont Change from here downward  ====================================

// Get Pixel size from the image
def cal = getCurrentServer().getPixelCalibration()
def pixelWidth = cal.pixelWidth
def pixelHeight = cal.pixelHeight
//print("============ pixelWidth="+ pixelWidth+ ", pixelHeight="+pixelHeight+" =====================")

var stardist = StarDist2D.builder(pathModel)
      .threshold(param_threshold)              // Prediction threshold
      .normalizePercentiles(normalize_low_pct,normalize_high_pct) // Percentile normalization
      .pixelSize(pixelWidth)              // Resolution for detection    
      .measureShape()
      .cellExpansion(param_expansion) //Cell expansion in microns
      .build()

// ===================  Segment Whole Tissue and Cells within it, convert Cells and Nuc from detections to annotationa  ====================================
var imageData = getCurrentImageData()
if (segmentTissue) {
    //createAnnotationsFromPixelClassifier(WholeTissueClassifier, WholeTissue_MinSize, WholeTissue_MinHoleSize, "SELECT_NEW")
    createAnnotationsFromPixelClassifier(WholeTissueClassifier, WholeTissue_MinSize, WholeTissue_MinHoleSize)
    
    println '======================== Whole Tissue segmentation Done =================== '
}
if (segmentCells) {

    resetSelection()
    selectObjectsByClassification(wholeTissueClass);

    var pathObjects = getSelectedObjects()
    if (pathObjects.isEmpty()) {
        Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
        return
    }
    stardist.detectObjects(imageData, pathObjects)

    // Get all the cells and turn them into annotations, so they can be edited
    def cellClass = getPathClass(cellClassName)
    //def nucClass = getPathClass(nucClassName)
    for (cell in getCellObjects()) {cell.setPathClass(cellClass) }
    
    println '======================== Cell segmentation Done =========================== '
}


// ===================  Load Visium HD Spots and associate them to Cells ====================================
if (loadSpots) {
    // Create BufferedReader
    def csvReader = new BufferedReader(new FileReader(csvfile));
    def plane = ImagePlane.getDefaultPlane();

    listOfObjects = [];    
    int first_row = 1;
    // Loop through all the rows of the CSV file.
    while ((row = csvReader.readLine()) != null) {
    
        if (first_row)
        {
            first_row = 0;
        }
        else 
        {
            def rowContent = row.split(",")
            String barcode = rowContent[0] as String;
            int  in_tissue = rowContent[1] as int;
            int  array_row = rowContent[2] as int;
            int  array_col = rowContent[3] as int;
            double cx = rowContent[4] as double;
            double cy = rowContent[5] as double;
        
            
            int first_time = 1;
            // Create annotation
            if (in_tissue) {
                double px = cx - spot_diameter_fullres/2;
                double py = cy - spot_diameter_fullres/2;
                def roi = new RectangleROI(py, px, spot_diameter_fullres, spot_diameter_fullres, plane)
                
                // Spots are imported as detections
                //def detection = new PathDetectionObject(roi, PathClass.fromString("Spot"));                
                def detection = new PathTileObject(roi, PathClass.fromString("Spot"),null);                
                if (first_time) {
                    first_time = 0;
                }
                detection.getMeasurementList().putMeasurement("array_row", array_row);
                detection.getMeasurementList().putMeasurement("array_col", array_col);
                detection.getMeasurementList().putMeasurement("cx", cx);
                detection.getMeasurementList().putMeasurement("cy", cy);
                detection.getMeasurementList().putMeasurement("InCell", 0);
                detection.getMeasurementList().putMeasurement("InNuc", 0);
                //detection.setName(barcode);
                detection.setName(barcode.replace('"',''));
                listOfObjects << detection;        
            }
        }
    }
    //imageData.getHierarchy().addObjects(listOfAnnotation, true);
    addObjects(listOfObjects);
    println '======================== Spots imported ==================================='
}

if (associateSpotsToCells) {
    def cellObj = getCellObjects() 
    for (def cell in cellObj) {
        def subCellObj = getCurrentHierarchy().getObjectsForROI(null, cell.getROI()) 
        nSpots = 0
        for (s in subCellObj) {
            if (s.isTile()) {
                nSpots++
                parentName = cell.getID().toString()
                sName = s.getName()
                //newName = sName + "__" + parentName;
                //newName = String.format("%s__%s", sName, parentName)
                newName = String.join("",sName, "__", parentName)
                s.setName(newName)
                s.measurements.put("InCell", 1)
                s.measurements.put("InNuc", 0)
            }        
        }
        def subNucObj = getCurrentHierarchy().getObjectsForROI(null, cell.getNucleusROI()) // get all objects within each cellObj Nucleus ROI
        nNucSpots = 0
        for (s in subNucObj) {
            if (s.isTile()) {
                nNucSpots++
                s.measurements.put("InNuc", 1)
            }
        }
        
        cell.measurements.put("nSpots",nSpots)
        cell.measurements.put("nNucSpots",nNucSpots)
    }
        
    println '======================== Associate Spots To Cells Done ==================== '
}


// ===================  Add pixel classifier measurements to spots ====================================
if (runPixelClassifierForSpot) {
    selectTiles()
    addPixelClassifierMeasurements(PixelClassifier, PixelClassifier)
    resetSelection()
}

// ===================  Add pixel classifier measurements to cells ====================================
if (runPixelClassifierForCell) {
    selectObjectsByClassification("Cell");
    addPixelClassifierMeasurements(PixelClassifier, PixelClassifier)
    resetSelection()
}

// ===================  Add Measurements to Cells  ====================================
if (AddMeasurementsToCells) {
    selectObjectsByClassification("Cell");
    runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons":2.0,"region":"ROI","tileSizeMicrons":25.0,"colorOD":false,"colorStain1":true,"colorStain2":true,"colorStain3":false,"colorRed":false,"colorGreen":false,"colorBlue":false,"colorHue":false,"colorSaturation":false,"colorBrightness":false,"doMean":true,"doStdDev":true,"doMinMax":true,"doMedian":true,"doHaralick":false,"haralickDistance":1,"haralickBins":32}')
    addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")
    selectAnnotations();
    runPlugin('qupath.lib.plugins.objects.SmoothFeaturesPlugin', '{"fwhmMicrons":25.0,"smoothWithinClasses":false}')
}

// ===================  Export Cells as GeoJson  ====================================
if (exportCellsAsGeoJson) {    
    
    File directory = new File(buildFilePath(PROJECT_BASE_DIR,resultsSubFolder));
    directory.mkdirs();
    imageName = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())

    cells = getCellObjects()
    selectObjects(cells)  
    exportSelectedObjectsToGeoJson(buildFilePath(directory.toString(),imageName+'_cells.geojson'), "EXCLUDE_MEASUREMENTS", "PRETTY_JSON", "FEATURE_COLLECTION")        
}


// ===================  Save Results ====================================
if (saveResultTable)
{        
    println '====================== Save Results Table ... ======================='
    File directory = new File(buildFilePath(PROJECT_BASE_DIR,resultsSubFolder));
    directory.mkdirs();
    imageName = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())
    saveAnnotationMeasurements(buildFilePath(directory.toString(),imageName+'_annotations.csv'));
    saveDetectionMeasurements(buildFilePath(directory.toString(),imageName+'_detections.csv'));
}

println '======================== Workflow Done! ==================================='