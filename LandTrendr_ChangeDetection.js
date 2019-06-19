//######################################################################################################## 
//#                                                                                                    #\\
//#                        LANDTRENDR GREATEST MAGNITUDE DISTURBANCE MAPPING                           #\\
//#                                                                                                    #\\
//########################################################################################################


// date: 2019-06-19
// author: Javier Lopatin | javier.lopatin@kit.edu


//########################################################################################################
//##### INPUTS ##### 
//########################################################################################################

// define a geometry - there are lots of ways to do this, see the GEE User guide
var aoi = table; // should be a GEE geometry object - here we are getting it from an drawn polygon


// define years and dates to include in landsat image collection
var startYear  = 1985;    // what year do you want to start the time series 
var endYear    = 2019;    // what year do you want to end the time series
var startDay   = '01-01'; // what is the beginning of date filter | month-day
var endDay     = '03-30'; // what is the end of date filter | month-day

// define disturbance mapping filter parameters 
var treeLoss1  = 175;      // delta filter for 1 year duration disturbance, <= will not be included as disturbance - units are in units of VI defined in the following function definition
var treeLoss20 = 200;      // delta filter for 20 year duration disturbance, <= will not be included as disturbance - units are in units of VI defined in the following function definition
var preVal     = 400;      // pre-disturbance value threshold - values below the provided threshold will exclude disturbance for those pixels - units are in units of VI defined in the following function definition
var mmu        = 10;       // minimum mapping unit for disturbance patches - units of pixels


// define function to calculate a spectral index to segment with LT
var NDVI = function(img) {
    var index = img.normalizedDifference(['B4', 'B3'])                      // calculate NDVI (B4-B3)/(B4+B3)
                   .multiply(1000)                                          // ...scale results by 1000 so we can convert to int and retain some precision
                   .select([0], ['NDVI'])                                   // ...name the band
                   .set('system:time_start', img.get('system:time_start')); // ...set the output system:time_start metadata to the input image time_start otherwise it is null
    return index ;
};

var GNDVI = function(img) {
    var index = img.normalizedDifference(['B4', 'B2'])  // calculate GNDVI (B4-B2)/(B4+B2)
                   .multiply(1000)                                          
                   .select([0], ['NDVI'])                                   
                   .set('system:time_start', img.get('system:time_start')); 
    return index ;
};

var LSWI = function(img) {
    var index = img.normalizedDifference(['B4', 'B7'])  // calculate land sourface water index (B4-B7)/(B4+B7)
                   .multiply(1000)                                          
                   .select([0], ['LSWI'])                                   
                   .set('system:time_start', img.get('system:time_start')); 
    return index ;
};

var EVI = function(img) {  
    var index = img.expression( // calculate EVI ( 2.5 * ((Band 4 – Band 3) / (Band 4 + 6 * Band 3 – 7.5 * Band 1 + 1)) )
      '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': img.select('B5'),
      'RED': img.select('B4'),
      'BLUE': img.select('B2')})                    
                   .multiply(1000)                                          
                   .select([0], ['EVI'])                                    
                   .set('system:time_start', img.get('system:time_start')); 
    return index ;
};

var NBR = function(img) {
    var index = img.normalizedDifference(['B4', 'B7']) // calculate normalized difference of band 4 and band 7 (B4-B7)/(B4+B7)
                   .multiply(1000)                                          
                   .select([0], ['NBR'])                                    
                   .set('system:time_start', img.get('system:time_start')); 
    return index ;
};

//########################################################################################################
// define the spectral index to use in the analisis
var VI = NBR

// define the sign of spectral delta for vegetation loss for the segmentation index - 
// VI delta is negetive for vegetation loss, so -1 for NDVI, EVI, GNDVI, etc
var distDir = -1; 

//########################################################################################################


// define the segmentation parameters:
// reference: Kennedy, R. E., Yang, Z., & Cohen, W. B. (2010). Detecting trends in forest disturbance and recovery using yearly Landsat time series: 1. LandTrendr—Temporal segmentation algorithms. Remote Sensing of Environment, 114(12), 2897-2910.
//            https://github.com/eMapR/LT-GEE
var run_params = { 
  maxSegments:            6,
  spikeThreshold:         0.9,
  vertexCountOvershoot:   3,
  preventOneYearRecovery: true,
  recoveryThreshold:      0.25,
  pvalThreshold:          0.05,
  bestModelProportion:    0.75,
  minObservationsNeeded:  6
};

//########################################################################################################
//##### ANNUAL SR TIME SERIES COLLECTION BUILDING FUNCTIONS ##### 
//########################################################################################################

//----- MAKE A DUMMY COLLECTOIN FOR FILLTING MISSING YEARS -----
var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]); // make an image collection from an image with 6 bands all set to 0 and then make them masked values


//------ L8 to L7 HARMONIZATION FUNCTION -----
// slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http://dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
var harmonizationRoy = function(oli) {
  var slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]);        // create an image of slopes per band for L8 TO L7 regression line - David Roy
  var itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]);     // create an image of y-intercepts per band for L8 TO L7 regression line - David Roy
  var y = oli.select(['B2','B3','B4','B5','B6','B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7']) // select OLI bands 2-7 and rename them to match L7 band names
             .resample('bicubic')                                                          // ...resample the L8 bands using bicubic
             .subtract(itcp.multiply(10000)).divide(slopes)                                // ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
             .set('system:time_start', oli.get('system:time_start'));                      // ...set the output system:time_start metadata to the input image time_start otherwise it is null
  return y.toShort();                                                                       // return the image as short to match the type of the other data
};


//------ RETRIEVE A SENSOR SR COLLECTION FUNCTION -----
var getSRcollection = function(year, startDay, endDay, sensor, aoi) {
  // get a landsat collection for given year, day range, and sensor
  var srCollection = ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR') // get surface reflectance images
                       .filterBounds(aoi)                                  // ...filter them by intersection with AOI
                       .filterDate(year+'-'+startDay, year+'-'+endDay);    // ...filter them by year and day range
  
  // apply the harmonization function to LC08 (if LC08), subset bands, unmask, and resample           
  srCollection = srCollection.map(function(img) {
    var dat = ee.Image(
      ee.Algorithms.If(
        sensor == 'LC08',                                                  // condition - if image is OLI
        harmonizationRoy(img.unmask()),                                    // true - then apply the L8 TO L7 alignment function after unmasking pixels that were previosuly masked (why/when are pixels masked)
        img.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])                   // false - else select out the reflectance bands from the non-OLI image
           .unmask()                                                       // ...unmask any previously masked pixels 
           .resample('bicubic')                                            // ...resample by bicubic 
           .set('system:time_start', img.get('system:time_start'))         // ...set the output system:time_start metadata to the input image time_start otherwise it is null
      )
    );
    
    // make a cloud, cloud shadow, and snow mask from fmask band
    var qa = img.select('pixel_qa');                                       // select out the fmask band
    var mask = qa.bitwiseAnd(8).eq(0).and(                                 // include shadow
               qa.bitwiseAnd(16).eq(0)).and(                               // include snow
               qa.bitwiseAnd(32).eq(0));                                   // include clouds
    
    // apply the mask to the image and return it
    return dat.mask(mask); //apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
  });

  return srCollection; // return the prepared collection
};


//------ FUNCTION TO COMBINE LT05, LE07, & LC08 COLLECTIONS -----
var getCombinedSRcollection = function(year, startDay, endDay, aoi) {
    var lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi);       // get TM collection for a given year, date range, and area
    var le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi);       // get ETM+ collection for a given year, date range, and area
    var lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi);       // get OLI collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8)); // merge the individual sensor collections into one imageCollection object
    return mergedCollection;                                              // return the Imagecollection
};


//------ FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE PER YEAR BY MEDOID -----
/*
  LT expects only a single image per year in a time series, there are lost of ways to
  do best available pixel compositing - we have found that a mediod composite requires little logic
  is robust, and fast
  
  Medoids are representative objects of a data set or a cluster with a data set whose average 
  dissimilarity to all the objects in the cluster is minimal. Medoids are similar in concept to 
  means or centroids, but medoids are always members of the data set.
*/

// make a medoid composite with equal weight among indices
var medoidMosaic = function(inCollection, dummyCollection) {
  
  // fill in missing years with the dummy collection
  var imageCount = inCollection.toList(1).length();                                                            // get the number of images 
  var finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection)); // if the number of images in this year is 0, then use the dummy collection, otherwise use the SR collection
  
  // calculate median across images in collection per band
  var median = finalCollection.median();                                                                       // calculate the median of the annual image collection - returns a single 6 band image - the collection median per band
  
  // calculate the different between the median and the observation per image per band
  var difFromMedian = finalCollection.map(function(img) {
    var diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));                                       // get the difference between each image/band and the corresponding band median and take to power of 2 to make negatives positive and make greater differences weight more
    return diff.reduce('sum').addBands(img);                                                                   // per image in collection, sum the powered difference across the bands - set this as the first band add the SR bands to it - now a 7 band image collection
  });
  
  // get the medoid by selecting the image pixel with the smallest difference between median and observation per band 
  return ee.ImageCollection(difFromMedian).reduce(ee.Reducer.min(7)).select([1,2,3,4,5,6], ['B1','B2','B3','B4','B5','B7']); // find the powered difference that is the least - what image object is the closest to the median of teh collection - and then subset the SR bands and name them - leave behind the powered difference band
};


//------ FUNCTION TO APPLY MEDOID COMPOSITING FUNCTION TO A COLLECTION -------------------------------------------
var buildMosaic = function(year, startDay, endDay, aoi, dummyCollection) {                                                                      // create a temp variable to hold the upcoming annual mosiac
  var collection = getCombinedSRcollection(year, startDay, endDay, aoi);  // get the SR collection
  var img = medoidMosaic(collection, dummyCollection)                     // apply the medoidMosaic function to reduce the collection to single image per year by medoid 
              .set('system:time_start', (new Date(year,8,1)).valueOf());  // add the year to each medoid image - the data is hard-coded Aug 1st 
  return ee.Image(img);                                                   // return as image object
};


//------ FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
var buildMosaicCollection = function(startYear, endYear, startDay, endDay, aoi, dummyCollection) {
  var imgs = [];                                                                    // create empty array to fill
  for (var i = startYear; i <= endYear; i++) {                                      // for each year from hard defined start to end build medoid composite and then add to empty img array
    var tmp = buildMosaic(i, startDay, endDay, aoi, dummyCollection);               // build the medoid mosaic for a given year
    imgs = imgs.concat(tmp.set('system:time_start', (new Date(i,8,1)).valueOf()));  // concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the year that is being worked on for Aug 1st
  }
  return ee.ImageCollection(imgs);                                                  // return the array img array as an image collection
};

//########################################################################################################
//##### UNPACKING LT-GEE OUTPUT STRUCTURE FUNCTIONS ##### 
//########################################################################################################

// ----- FUNCTION TO EXTRACT VERTICES FROM LT RESULTS AND STACK BANDS -----
var getLTvertStack = function(LTresult) {
  var emptyArray = [];                              // make empty array to hold another array whose length will vary depending on maxSegments parameter    
  var vertLabels = [];                              // make empty array to hold band names whose length will vary depending on maxSegments parameter 
  var iString;                                      // initialize variable to hold vertex number
  for(var i=1;i<=run_params.maxSegments+1;i++){     // loop through the maximum number of vertices in segmentation and fill empty arrays
    iString = i.toString();                         // define vertex number as string 
    vertLabels.push("vert_"+iString);               // make a band name for given vertex
    emptyArray.push(0);                             // fill in emptyArray
  }
  
  var zeros = ee.Image(ee.Array([emptyArray,        // make an image to fill holes in result 'LandTrendr' array where vertices found is not equal to maxSegments parameter plus 1
                                 emptyArray,
                                 emptyArray]));
  
  var lbls = [['yrs_','src_','fit_'], vertLabels,]; // labels for 2 dimensions of the array that will be cast to each other in the final step of creating the vertice output 

  var vmask = LTresult.arraySlice(0,3,4);           // slices out the 4th row of a 4 row x N col (N = number of years in annual stack) matrix, which identifies vertices - contains only 0s and 1s, where 1 is a vertex (referring to spectral-temporal segmentation) year and 0 is not
  
  var ltVertStack = LTresult.arrayMask(vmask)       // uses the sliced out isVert row as a mask to only include vertice in this data - after this a pixel will only contain as many "bands" are there are vertices for that pixel - min of 2 to max of 7. 
                      .arraySlice(0, 0, 3)          // ...from the vertOnly data subset slice out the vert year row, raw spectral row, and fitted spectral row
                      .addBands(zeros)              // ...adds the 3 row x 7 col 'zeros' matrix as a band to the vertOnly array - this is an intermediate step to the goal of filling in the vertOnly data so that there are 7 vertice slots represented in the data - right now there is a mix of lengths from 2 to 7
                      .toArray(1)                   // ...concatenates the 3 row x 7 col 'zeros' matrix band to the vertOnly data so that there are at least 7 vertice slots represented - in most cases there are now > 7 slots filled but those will be truncated in the next step
                      .arraySlice(1, 0, run_params.maxSegments+1) // ...before this line runs the array has 3 rows and between 9 and 14 cols depending on how many vertices were found during segmentation for a given pixel. this step truncates the cols at 7 (the max verts allowed) so we are left with a 3 row X 7 col array
                      .arrayFlatten(lbls, '');      // ...this takes the 2-d array and makes it 1-d by stacking the unique sets of rows and cols into bands. there will be 7 bands (vertices) for vertYear, followed by 7 bands (vertices) for rawVert, followed by 7 bands (vertices) for fittedVert, according to the 'lbls' list

  return ltVertStack;                               // return the stack
};


//########################################################################################################
//##### GREATEST DISTURBANCE EXTRACTION FUNCTIONS #####
//########################################################################################################

// ----- function to extract greatest disturbance based on spectral delta between vertices 
var extractDisturbance = function(lt, distDir, params, mmu) {
  // select only the vertices that represents a change
  var vertexMask = lt.arraySlice(0, 3, 4); // get the vertex - yes(1)/no(0) dimension
  var vertices = lt.arrayMask(vertexMask); // convert the 0's to masked
  
  // construct segment start and end point years and index values
  var left = vertices.arraySlice(1, 0, -1);    // slice out the vertices as the start of segments
  var right = vertices.arraySlice(1, 1, null); // slice out the vertices as the end of segments
  var startYear = left.arraySlice(0, 0, 1);    // get year dimension of LT data from the segment start vertices
  var startVal = left.arraySlice(0, 2, 3);     // get spectral index dimension of LT data from the segment start vertices
  var endYear = right.arraySlice(0, 0, 1);     // get year dimension of LT data from the segment end vertices 
  var endVal = right.arraySlice(0, 2, 3);      // get spectral index dimension of LT data from the segment end vertices
  
  var dur = endYear.subtract(startYear);       // subtract the segment start year from the segment end year to calculate the duration of segments 
  var mag = endVal.subtract(startVal);         // substract the segment start index value from the segment end index value to calculate the delta of segments 

  // concatenate segment start year, delta, duration, and starting spectral index value to an array 
  var distImg = ee.Image.cat([startYear.add(1), mag, dur, startVal.multiply(distDir)]).toArray(0); // make an image of segment attributes - multiply by the distDir parameter to re-orient the spectral index if it was flipped for segmentation - do it here so that the subtraction to calculate segment delta in the above line is consistent - add 1 to the detection year, because the vertex year is not the first year that change is detected, it is the following year
 
  // sort the segments in the disturbance attribute image delta by spectral index change delta  
  var distImgSorted = distImg.arraySort(mag.multiply(-1));                                  // flip the delta around so that the greatest delta segment is first in order

  // slice out the first (greatest) delta
  var tempDistImg = distImgSorted.arraySlice(1, 0, 1).unmask(ee.Image(ee.Array([[0],[0],[0],[0]])));           // get the first segment in the sorted array

  // make an image from the array of attributes for the greatest disturbance
  var finalDistImg = ee.Image.cat(tempDistImg.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['yod']]),     // slice out year of disturbance detection and re-arrange to an image band 
                                  tempDistImg.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['mag']]),     // slice out the disturbance magnitude and re-arrange to an image band 
                                  tempDistImg.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['dur']]),     // slice out the disturbance duration and re-arrange to an image band
                                  tempDistImg.arraySlice(0,3,4).arrayProject([1]).arrayFlatten([['preval']])); // slice out the pre-disturbance spectral value and re-arrange to an image band
  
  // filter out disturbances based on user settings
  var threshold = ee.Image(finalDistImg.select(['dur']))                        // get the disturbance band out to apply duration dynamic disturbance magnitude threshold 
                    .multiply((params.tree_loss20 - params.tree_loss1) / 19.0)  // ...
                    .add(params.tree_loss1)                                     //    ...interpolate the magnitude threshold over years between a 1-year mag thresh and a 20-year mag thresh
                    .lte(finalDistImg.select(['mag']))                          // ...is disturbance less then equal to the interpolated, duration dynamic disturbance magnitude threshold 
                    .and(finalDistImg.select(['mag']).gt(0))                    // and is greater than 0  
                    .and(finalDistImg.select(['preval']).gt(params.pre_val));   // and is greater than pre-disturbance spectral index value threshold
  
  // apply the filter mask
  finalDistImg = finalDistImg.mask(threshold).int16(); 
  
   // patchify the remaining disturbance pixels using a minimum mapping unit
  if(mmu > 1){
    var mmuPatches = finalDistImg.select(['yod'])           // patchify based on disturbances having the same year of detection
                            .connectedPixelCount(mmu, true) // count the number of pixel in a candidate patch
                            .gte(mmu);                      // are the the number of pixels per candidate patch greater than user-defined minimum mapping unit?
    finalDistImg = finalDistImg.updateMask(mmuPatches);     // mask the pixels/patches that are less than minimum mapping unit
  } 
  
  return finalDistImg; // return the filtered greatest disturbance attribute image
};


//########################################################################################################
//##### BUILD COLLECTION AND RUN LANDTRENDR #####
//########################################################################################################

//----- BUILD LT COLLECTION -----
// build annual surface reflection collection
var annualSRcollection = buildMosaicCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection); // put together the cloud-free medoid surface reflectance annual time series collection

// apply the function to calculate the segmentation index and adjust the values by the distDir parameter - flip index so that a vegetation loss is associated with a postive delta in spectral value
var ltCollection = annualSRcollection.map(VI)                                             // map the function over every image in the collection - returns a 1-band annual image collection of the spectral index
                                              .map(function(img) {return img.multiply(distDir)           // ...multiply the segmentation index by the distDir to ensure that vegetation loss is associated with a positive spectral delta
                                              .set('system:time_start', img.get('system:time_start'))}); // ...set the output system:time_start metadata to the input image time_start otherwise it is null

//----- RUN LANDTRENDR -----
run_params.timeSeries = ltCollection;               // add LT collection to the segmentation run parameter object
var lt = ee.Algorithms.TemporalSegmentation.LandTrendr(run_params); // run LandTrendr spectral temporal segmentation algorithm



//########################################################################################################
//##### RUN THE GREATEST DISTURBANCE EXTRACT FUCTION #####
//########################################################################################################

// assemble the disturbance extraction parameters
var distParams = {
  tree_loss1: treeLoss1,
  tree_loss20: treeLoss20,  
  pre_val: preVal           
};

// run the dist extract function
var distImg = extractDisturbance(lt.select('LandTrendr'), distDir, distParams).clip(aoi);


//########################################################################################################
//##### DISTURBANCE MAP DISPLAY #####
//########################################################################################################

// ----- set visualization dictionaries -----
var yodVizParms = {
  min: startYear+1,
  max: endYear,
  palette: ['#9400D3', '#4B0082', '#0000FF', '#00FF00', '#FFFF00', '#FF7F00', '#FF0000']
};

var magVizParms = {
  min: distParams.tree_loss1,
  max: 1000,
  palette: ['#0000FF', '#00FF00', '#FFFF00', '#FF7F00', '#FF0000']
};

var durVizParms = {
  min: 1,
  max: endYear-startYear,
  palette: ['#FF0000', '#FF7F00', '#FFFF00', '#00FF00', '#0000FF']
};

var preValVizParms = {
  min: preVal,
  max: 800,
  palette: ['#FF0000', '#FF7F00', '#FFFF00', '#00FF00', '#0000FF']
};


// ----- display the disturbance attribute maps ----- 
Map.centerObject(aoi, 12);                                                    // center the map display and set zoom level
//Map.addLayer(distImg.select(['preval']), preValVizParms, 'Pre-dist Value'); // add pre-disturbacne spectral index value to map
//Map.addLayer(distImg.select(['dur']), durVizParms, 'Duration');             // add disturbance duration to map
//Map.addLayer(distImg.select(['mag']), magVizParms, 'Magnitude');            // add magnitude to map
Map.addLayer(distImg.select(['yod']), yodVizParms, 'Year of Detection');      // add disturbance year of detection to map

// Save results to Drive
Export.image.toDrive({
  image: distImg.select(['preval']),
  description: 'LandTrendr_preval',
  folder: 'data_earth_engine',
  scale: 30,
  region: aoi
});
Export.image.toDrive({
  image: distImg.select(['dur']),
  description: 'LandTrendr_duration',
  folder: 'data_earth_engine',
  scale: 30,
  region: aoi
});
Export.image.toDrive({
  image: distImg.select(['mag']),
  description: 'LandTrendr_magnitude',
  folder: 'data_earth_engine',
  scale: 30,
  region: aoi
});
Export.image.toDrive({
  image: distImg.select(['yod']),
  description: 'LandTrendr_year-of-detection',
  folder: 'data_earth_engine',
  scale: 30,
  region: aoi
});
