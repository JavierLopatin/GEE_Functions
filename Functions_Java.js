//########################################################################################################
//#                                                                                                    #\\
//#                                         Set of GEE Java-based funcitons                                        #\\
//#                                                                                                    #\\
//########################################################################################################


// date: 2018-11
// author: Javier Lopatin | javier.lopatin@kit.edu; javierlopatin@gmail.com

// With contributions of: - https://github.com/eMapR/LT-GEE



//########################################################################################################
//##### ANNUAL SR TIME SERIES COLLECTION BUILDING FUNCTIONS #####
//########################################################################################################


//------ L8 to L7 HARMONIZATION FUNCTION -----
// slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F.,
// Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8
// reflective wavelength and normalized difference vegetation index continuity,
// Remote Sensing of Environment, 185, 57-70.(http://dx.doi.org/10.1016/j.rse.2015.12.024);
// Table 2 - reduced major axis (RMA) regression coefficients
var harmonizationRoy = function(oli) {
  var slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]);        // RMA - create an image of slopes per band for L8 TO L7 regression line - David Roy
  var itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]);     // RMA - create an image of y-intercepts per band for L8 TO L7 regression line - David Roy
  //var slopes = ee.Image.constant([0.885, 0.9317, 0.9372, 0.8339, 0.8639, 0.9165]);  // least squares OLI to ETM+
  //var itcp = ee.Image.constant([0.0183, 0.0123, 0.0123, 0.0448, 0.0306, 0.0116]);   // least squares OLI to ETM+
  var y = oli.select(['B2','B3','B4','B5','B6','B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7']) // select OLI bands 2-7 and rename them to match L7 band names
             .resample('bicubic')                                                          // ...resample the L8 bands using bicubic
             .subtract(itcp.multiply(10000)).divide(slopes)                                // ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
             .set('system:time_start', oli.get('system:time_start'));                      // ...set the output system:time_start metadata to the input image time_start otherwise it is null
  return y.toShort();                                                                       // return the image as short to match the type of the other data
};
exports.harmonizationRoy = harmonizationRoy;


//------ FILTER A COLLECTION FUNCTION -----
var filterCollection = function(year, startDay, endDay, sensor, aoi){
  return ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR')
           .filterBounds(aoi)
           .filterDate(year+'-'+startDay, year+'-'+endDay);
};


//------ BUILD A COLLECTION FOR A GIVEN SENSOR AND YEAR -----
var buildSensorYearCollection = function(year, startDay, endDay, sensor, aoi){
  var startMonth = parseInt(startDay.substring(0, 2));
  var endMonth = parseInt(endDay.substring(0, 2));
  var srCollection;
  if(startMonth > endMonth){
    var oldYear = (parseInt(year)-1).toString();
    var newYear = year;
    var oldYearStartDay = startDay;
    var oldYearEndDay = '12-31';
    var newYearStartDay = '01-01';
    var newYearEndDay = endDay;

    var oldYearCollection = filterCollection(oldYear, oldYearStartDay, oldYearEndDay, sensor, aoi);
    var newYearCollection = filterCollection(newYear, newYearStartDay, newYearEndDay, sensor, aoi);

    srCollection = ee.ImageCollection(oldYearCollection.merge(newYearCollection));
  } else {
    srCollection = filterCollection(year, startDay, endDay, sensor, aoi);
  }

  return srCollection;
};
exports.buildSensorYearCollection = buildSensorYearCollection



//------ RETRIEVE A SENSOR SR COLLECTION FUNCTION -----
var getSRcollection = function(year, startDay, endDay, sensor, aoi, maskThese) {
  // make sure that mask labels are correct
  maskThese = (typeof maskThese !== 'undefined') ?  maskThese : ['cloud','shadow','snow'];
  var maskOptions = ['cloud', 'shadow', 'snow', 'water'];
  for(var i in maskThese){
    maskThese[i] =  maskThese[i].toLowerCase();
    var test = maskOptions.indexOf(maskThese[i]);
    if(test == -1){
      print('error: '+maskThese[i]+' is not included in the list maskable features. Please see ___ for list of maskable features to include in the maskThese parameter');
      return 'error';
    }
  }

  // get a landsat collection for given year, day range, and sensor
  var srCollection = buildSensorYearCollection(year, startDay, endDay, sensor, aoi);

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

//    // make a cloud, cloud shadow, and snow mask from fmask band
//    var qa = img.select('pixel_qa');                                       // select out the fmask band
//    var mask = qa.bitwiseAnd(8).eq(0).and(                                 // include shadow
//               qa.bitwiseAnd(16).eq(0)).and(                               // include snow
//               qa.bitwiseAnd(32).eq(0));                                   // include clouds

    var mask = ee.Image(1);
    if(maskThese.length !== 0){
      var qa = img.select('pixel_qa');
      for(var i in maskThese){
        if(maskThese[i] == 'water'){mask = qa.bitwiseAnd(4).eq(0).multiply(mask)}
        if(maskThese[i] == 'shadow'){mask = qa.bitwiseAnd(8).eq(0).multiply(mask)}
        if(maskThese[i] == 'snow'){mask = qa.bitwiseAnd(16).eq(0).multiply(mask)}
        if(maskThese[i] == 'cloud'){mask = qa.bitwiseAnd(32).eq(0).multiply(mask)}
      }
      return dat.mask(mask); //apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
    } else{
      return dat;
    }
  });

  return srCollection; // return the prepared collection
};
exports.getSRcollection = getSRcollection;


//------ FUNCTION TO COMBINE LT05, LE07, & LC08 COLLECTIONS -----
var getCombinedSRcollection = function(year, startDay, endDay, aoi, maskThese) {
    var lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi, maskThese);       // get TM collection for a given year, date range, and area
    var le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi, maskThese);       // get ETM+ collection for a given year, date range, and area
    var lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi, maskThese);       // get OLI collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8)); // merge the individual sensor collections into one imageCollection object
    return mergedCollection;                                              // return the Imagecollection
};
exports.getCombinedSRcollection = getCombinedSRcollection



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
var buildMosaic = function(year, startDay, endDay, aoi, dummyCollection, maskThese) {                                                                      // create a temp variable to hold the upcoming annual mosiac
  var collection = getCombinedSRcollection(year, startDay, endDay, aoi, maskThese);  // get the SR collection
  var img = medoidMosaic(collection, dummyCollection)                     // apply the medoidMosaic function to reduce the collection to single image per year by medoid
              .set('system:time_start', (new Date(year,8,1)).valueOf());  // add the year to each medoid image - the data is hard-coded Aug 1st
  return ee.Image(img);                                                   // return as image object
};


//------ FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
var buildSRcollection = function(startYear, endYear, startDay, endDay, aoi, maskThese) {
  var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]); // make an image collection from an image with 6 bands all set to 0 and then make them masked values
  var imgs = [];                                                                         // create empty array to fill
  for (var i = startYear; i <= endYear; i++) {                                           // for each year from hard defined start to end build medoid composite and then add to empty img array
    var tmp = buildMosaic(i, startDay, endDay, aoi, dummyCollection, maskThese);                    // build the medoid mosaic for a given year
    imgs = imgs.concat(tmp.set('system:time_start', (new Date(i,8,1)).valueOf()));       // concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the year that is being worked on for Aug 1st
  }
  return ee.ImageCollection(imgs);                                                       // return the array img array as an image collection
};
exports.buildSRcollection = buildSRcollection;

//------ FUNCTION TO RETURN A LIST OF IMAGES THAT GO INTO ANNUAL SR COMPOSITE COLLECTION ------------------------------
function getImgIDs(img){return ee.String(ee.Image(img).get('system:id'));}

var getCollectionIDlist = function(startYear, endYear, startDay, endDay, aoi) {
  var first = true;
  for (var i = startYear; i <= endYear; i++){
    var lt5 = buildSensorYearCollection(i, startDay, endDay, 'LT05', aoi);
    var le7 = buildSensorYearCollection(i, startDay, endDay, 'LE07', aoi);
    var lc8 = buildSensorYearCollection(i, startDay, endDay, 'LC08', aoi);
    var tmp = ee.ImageCollection(lt5.merge(le7).merge(lc8));
    if(first === true){
      var all = tmp;
      first = false;
    } else{
      all = all.merge(tmp);
    }
  }
  return all.toList(all.size().add(1)).map(getImgIDs);
};
exports.getCollectionIDlist = getCollectionIDlist;



//------ FUNCTION TO COUNT NUMBER OF UNMASKED PIXELS IN AN INTRA ANNUAL COLLECTION ------------------------------
var countClearViewPixels = function(intraAnnualSRcollection){
  var binary = intraAnnualSRcollection.map(function(img){
    return img.select(0)
              .multiply(0)
              .add(1)
              .unmask(0);
  });
  return binary.sum();
};
exports.countClearViewPixels = countClearViewPixels;

//------ FUNCTION TO BUILD ANNAUL COLLECTION OF NUMBER OF UNMASKED PIXELS AVAILABLE TO BUILD COMPOSITE ------------------------------
var buildClearPixelCountCollection = function(startYear, endYear, startDay, endDay, aoi, maskThese) {
  var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]);
  var imgs = [];
  for (var i = startYear; i <= endYear; i++) {
    var collection = getCombinedSRcollection(i, startDay, endDay, aoi, maskThese, maskThese);
    var imageCount = collection.toList(1).length();
    var finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), collection, dummyCollection));
    var notMaskCount = countClearViewPixels(finalCollection);
    imgs = imgs.concat(notMaskCount.set('system:time_start', (new Date(i,8,1)).valueOf()));
  }
  return ee.ImageCollection(imgs);
};
exports.buildClearPixelCountCollection = buildClearPixelCountCollection


//########################################################################################################
//##### UNPACKING LT-GEE OUTPUT STRUCTURE FUNCTIONS #####
//########################################################################################################

// ----- FUNCTION TO EXTRACT VERTICES FROM LT RESULTS AND STACK BANDS -----
var getLTvertStack = function(lt, runParams) {
  lt = lt.select('LandTrendr');
  var emptyArray = [];                              // make empty array to hold another array whose length will vary depending on maxSegments parameter
  var vertLabels = [];                              // make empty array to hold band names whose length will vary depending on maxSegments parameter
  for(var i=1;i<=runParams.maxSegments+1;i++){     // loop through the maximum number of vertices in segmentation and fill empty arrays                        // define vertex number as string
    vertLabels.push("vert_"+i.toString());               // make a band name for given vertex
    emptyArray.push(0);                             // fill in emptyArray
  }

  var zeros = ee.Image(ee.Array([emptyArray,        // make an image to fill holes in result 'LandTrendr' array where vertices found is not equal to maxSegments parameter plus 1
                                 emptyArray,
                                 emptyArray]));

  var lbls = [['yrs_','src_','fit_'], vertLabels,]; // labels for 2 dimensions of the array that will be cast to each other in the final step of creating the vertice output

  var vmask = lt.arraySlice(0,3,4);           // slices out the 4th row of a 4 row x N col (N = number of years in annual stack) matrix, which identifies vertices - contains only 0s and 1s, where 1 is a vertex (referring to spectral-temporal segmentation) year and 0 is not

  var ltVertStack = lt.arrayMask(vmask)       // uses the sliced out isVert row as a mask to only include vertice in this data - after this a pixel will only contain as many "bands" are there are vertices for that pixel - min of 2 to max of 7.
                      .arraySlice(0, 0, 3)          // ...from the vertOnly data subset slice out the vert year row, raw spectral row, and fitted spectral row
                      .addBands(zeros)              // ...adds the 3 row x 7 col 'zeros' matrix as a band to the vertOnly array - this is an intermediate step to the goal of filling in the vertOnly data so that there are 7 vertice slots represented in the data - right now there is a mix of lengths from 2 to 7
                      .toArray(1)                   // ...concatenates the 3 row x 7 col 'zeros' matrix band to the vertOnly data so that there are at least 7 vertice slots represented - in most cases there are now > 7 slots filled but those will be truncated in the next step
                      .arraySlice(1, 0, runParams.maxSegments+1) // ...before this line runs the array has 3 rows and between 9 and 14 cols depending on how many vertices were found during segmentation for a given pixel. this step truncates the cols at 7 (the max verts allowed) so we are left with a 3 row X 7 col array
                      .arrayFlatten(lbls, '');      // ...this takes the 2-d array and makes it 1-d by stacking the unique sets of rows and cols into bands. there will be 7 bands (vertices) for vertYear, followed by 7 bands (vertices) for rawVert, followed by 7 bands (vertices) for fittedVert, according to the 'lbls' list

  return ltVertStack;                               // return the stack
};

exports.getLTvertStack = getLTvertStack;










// #######################################################################################
// ###### INDEX CALCULATION FUNCTIONS ####################################################
// #######################################################################################

// TASSELLED CAP
var tcTransform = function(img){
  var b = ee.Image(img).select(["B1", "B2", "B3", "B4", "B5", "B7"]); // select the image bands
  var brt_coeffs = ee.Image.constant([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]); // set brt coeffs - make an image object from a list of values - each of list element represents a band
  var grn_coeffs = ee.Image.constant([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]); // set grn coeffs - make an image object from a list of values - each of list element represents a band
  var wet_coeffs = ee.Image.constant([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]); // set wet coeffs - make an image object from a list of values - each of list element represents a band

  var sum = ee.Reducer.sum(); // create a sum reducer to be applyed in the next steps of summing the TC-coef-weighted bands
  var brightness = b.multiply(brt_coeffs).reduce(sum); // multiply the image bands by the brt coef and then sum the bands
  var greenness = b.multiply(grn_coeffs).reduce(sum); // multiply the image bands by the grn coef and then sum the bands
  var wetness = b.multiply(wet_coeffs).reduce(sum); // multiply the image bands by the wet coef and then sum the bands
  var angle = (greenness.divide(brightness)).atan().multiply(180/Math.PI).multiply(100);
  var tc = brightness.addBands(greenness)
                     .addBands(wetness)
                     .addBands(angle)
                     .select([0,1,2,3], ['TCB','TCG','TCW','TCA']) //stack TCG and TCW behind TCB with .addBands, use select() to name the bands
                     .set('system:time_start', img.get('system:time_start'));
  return tc;
};

// NBR
var nbrTransform = function(img) {
    var nbr = img.normalizedDifference(['B4', 'B7']) // calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
                 .multiply(1000) // scale results by 1000
                 .select([0], ['NBR']) // name the band
                 .set('system:time_start', img.get('system:time_start'));
    return nbr;
};

// NDVI
var ndviTransform = function(img){
  var ndvi = img.normalizedDifference(['B4', 'B3']) // calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
                .multiply(1000) // scale results by 1000
                .select([0], ['NDVI']) // name the band
                .set('system:time_start', img.get('system:time_start'));
  return ndvi;
};
                
// NDSI
var ndsiTransform = function(img){
  var ndsi = img.normalizedDifference(['B2', 'B5']) // calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
                .multiply(1000) // scale results by 1000
                .select([0], ['NDSI']) // name the band
                .set('system:time_start', img.get('system:time_start'));
  return ndsi;
};

// NDMI
var ndmiTransform = function(img) {
    var ndmi = img.normalizedDifference(['B4', 'B5']) // calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
                 .multiply(1000) // scale results by 1000
                 .select([0], ['NDMI']) // name the band
                 .set('system:time_start', img.get('system:time_start'));
    return ndmi;
};
// EVI
var eviTrabsform = funciton(img){
  var evi = img.expression(
      '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': img.select('B5'),
      'RED': img.select('B4'),
      'BLUE': img.select('B2')})
                   .multiply(1000)
                   .select([0], ['EVI'])
                   .set('system:time_start', img.get('system:time_start'));
    return evi ;
}

// CALCULATE A GIVEN INDEX
var calcIndex = function(img, index, flip){
  // make sure index string in upper case
  index = index.toUpperCase();

  // figure out if we need to calc tc
  var tcList = ['TCB', 'TCG', 'TCW', 'TCA'];
  var doTC = tcList.indexOf(index);
  if(doTC >= 0){
    var tc = tcTransform(img);
  }

  // need to flip some indices if this is intended for segmentation
  var indexFlip = 1;
  if(flip == 1){
    indexFlip = -1;
  }

  // need to cast raw bands to float to make sure that we don't get errors regarding incompatible bands
  // ...derivations are already float because of division or multiplying by decimal
  var indexImg;
  switch (index){
    case 'B1':
      indexImg = img.select(['B1']).float();//.multiply(indexFlip);
      break;
    case 'B2':
      indexImg = img.select(['B2']).float();//.multiply(indexFlip);
      break;
    case 'B3':
      indexImg = img.select(['B3']).float();//.multiply(indexFlip);
      break;
    case 'B4':
      indexImg = img.select(['B4']).multiply(indexFlip).float();
      break;
    case 'B5':
      indexImg = img.select(['B5']).float();//.multiply(indexFlip);
      break;
    case 'B7':
      indexImg = img.select(['B7']).float();//.multiply(indexFlip);
      break;
    case 'NBR':
      indexImg = nbrTransform(img).multiply(indexFlip);
      break;
    case 'NDMI':
      indexImg = ndmiTransform(img).multiply(indexFlip);
      break;
    case 'NDVI':
      indexImg = ndviTransform(img).multiply(indexFlip);
      break;
      case 'EVI':
        indexImg = eviTransform(img).multiply(indexFlip);
        break;
    case 'NDSI':
      indexImg = ndsiTransform(img).multiply(indexFlip);
      break;
    case 'TCB':
      indexImg = tc.select(['TCB'])//.multiply(indexFlip);
      break;
    case 'TCG':
      indexImg = tc.select(['TCG']).multiply(indexFlip);
      break;
    case 'TCW':
      indexImg = tc.select(['TCW']).multiply(indexFlip);
      break;
    case 'TCA':
      indexImg = tc.select(['TCA']).multiply(indexFlip);
      break;
    default:
      print('The index you provided is not supported');
  }

  return indexImg.set('system:time_start', img.get('system:time_start'));
};


exports.calcIndex = calcIndex;

// MAKE AN LT STACK
//var makeLTstack = function(img){
//  var allStack = calcIndex(img, index, 1); // calc index for segmentation
//  var ftvimg;
//  for(var ftv in ftvList){
//    ftvimg = calcIndex(img, ftvList[ftv], 0) // calc index for FTV
//                      .select([ftvList[ftv]],['ftv_'+ftvList[ftv].toLowerCase()]);
//
//    allStack = allStack.addBands(ftvimg)
//                       .set('system:time_start', img.get('system:time_start'));
//  }
//
//  return allStack;
//};

// arrange collection as an annual stack
// TODO: would be nice to name the bands with the original band name and the year
var TScollectionToStack = function(collection, startYear, endYear){
  var collectionArray = collection.toArrayPerBand();
  var nBands = ee.Image(collection.first()).bandNames().getInfo().length;//;
  var bandNames = getYearBandNames(startYear, endYear);
  var allStack = ee.Image();
  for (var i = 0; i < nBands; i++){
    var bandTS = collectionArray.select([i]).arrayFlatten([bandNames]);
    allStack = ee.Image.cat([allStack, bandTS]);
  }

  return allStack.slice(1,null);
};

exports.TScollectionToStack = TScollectionToStack;

var standardize = function(collection){
  var mean = collection.reduce(ee.Reducer.mean());
  var stdDev = collection.reduce(ee.Reducer.stdDev());

  var meanAdj = collection.map(function(img){
    return img.subtract(mean).set('system:time_start', img.get('system:time_start'));
  });

  return meanAdj.map(function(img){
    return img.divide(stdDev).set('system:time_start', img.get('system:time_start'));
  });
};


// STANDARDIZE TASSELED CAP BRIGHTNESS GREENNESS WETNESS AND REDUCE THE COLLECTION
var makeTCcomposite = function(annualSRcollection, reducer){
  var TCcomposite = annualSRcollection.map(function(img){
    var tcb = calcIndex(img, 'TCB', 1);//.unmask(0);
    var tcg = calcIndex(img, 'TCG', 1);//.unmask(0);
    var tcw = calcIndex(img, 'TCW', 1);//.unmask(0);
    return tcb.addBands(tcg)
              .addBands(tcw)
              .set('system:time_start', img.get('system:time_start'));
  });

  var tcb = TCcomposite.select(['TCB']);
  var tcg = TCcomposite.select(['TCG']);
  var tcw = TCcomposite.select(['TCW']);

  // standardize the TC bands
  var tcbStandard = standardize(tcb);
  var tcgStandard = standardize(tcg);
  var tcwStandard = standardize(tcw);

  // combine the standardized TC band collections into a single collection
  var tcStandard = tcbStandard.combine(tcgStandard).combine(tcwStandard);


  TCcomposite = tcStandard.map(function(img){
    var imgCollection = ee.ImageCollection.fromImages(
      [
        img.select(['TCB'],['Z']),
        img.select(['TCG'],['Z']),
        img.select(['TCW'],['Z'])
      ]
    );
    var reducedImg;
    switch(reducer){
      case 'mean':
        reducedImg = imgCollection.mean();
        break;
      case 'max':
        reducedImg = imgCollection.max();
        break;
      case 'sum':
        reducedImg = imgCollection.sum();
        break;
      default:
        print('The reducer you provided is not supported');
    }

    return reducedImg.multiply(1000).set('system:time_start', img.get('system:time_start'));
  });


  return TCcomposite;
};

// STANDARDIZE B5, TCB, TCG, NBR AND REDUCE THE COLLECTION
var makeEnsemblecomposite = function(annualSRcollection, reducer){
  // make a collection of the ensemble indices stacked as bands
  var stack = annualSRcollection.map(function(img){
    var b5   = calcIndex(img, 'B5', 1);
    var b7   = calcIndex(img, 'B7', 1);
    var tcw  = calcIndex(img, 'TCW', 1);
    var tca  = calcIndex(img, 'TCA', 1);
    var ndmi = calcIndex(img, 'NDMI', 1);
    var nbr  = calcIndex(img, 'NBR', 1);
    return b5.addBands(b7)
             .addBands(tcw)
             .addBands(tca)
             .addBands(ndmi)
             .addBands(nbr)
             .set('system:time_start', img.get('system:time_start'));
  });

  // make subset collections of each index
  var b5 = stack.select('B5');
  var b7 = stack.select('B7');
  var tcw = stack.select('TCW');
  var tca = stack.select('TCA');
  var ndmi = stack.select('NDMI');
  var nbr = stack.select('NBR');

  // standardize each index to mean 0 stdev 1
  var b5Standard = standardize(b5);
  var b7Standard = standardize(b7);
  var tcwStandard = standardize(tcw);
  var tcaStandard = standardize(tca);
  var ndmiStandard = standardize(ndmi);
  var nbrStandard = standardize(nbr);

  // combine the standardized band collections into a single collection
  var standard = b5Standard.combine(b7Standard).combine(tcwStandard).combine(tcaStandard)
                             .combine(ndmiStandard).combine(nbrStandard);

  // reduce the collection to a single value
  var composite = standard.map(function(img){
    var imgCollection = ee.ImageCollection.fromImages(
      [
        img.select(['B5'],['Z']),
        img.select(['B7'],['Z']),
        img.select(['TCW'],['Z']),
        img.select(['TCA'],['Z']),
        img.select(['NDMI'],['Z']),
        img.select(['NBR'],['Z']),
      ]
    );
    var reducedImg;
    switch(reducer){
      case 'mean':
        reducedImg = imgCollection.mean();
        break;
      case 'max':
        reducedImg = imgCollection.max();
        break;
      case 'sum':
        reducedImg = imgCollection.sum();
        break;
      default:
        print('The reducer you provided is not supported');
    }

    return reducedImg.multiply(1000).set('system:time_start', img.get('system:time_start'));
  });

  return composite;
};


makeEnsemblecomposite.exports = makeEnsemblecomposite;

// STANDARDIZE B5, TCB, TCG, NBR AND REDUCE THE COLLECTION
var makeEnsemblecomposite1 = function(annualSRcollection, reducer){
  // make a collection of the ensemble indices stacked as bands
  var TCcomposite = annualSRcollection.map(function(img){
    var b5   = calcIndex(img, 'B5', 1);
    var tcb  = calcIndex(img, 'TCB', 1);
    var tcg  = calcIndex(img, 'TCG', 1);
    var nbr  = calcIndex(img, 'NBR', 1);
    return b5.addBands(tcb)
             .addBands(tcg)
             .addBands(nbr)
             .set('system:time_start', img.get('system:time_start'));
  });

  // make subset collections of each index
  var b5 = TCcomposite.select('B5');
  var tcb = TCcomposite.select('TCB');
  var tcg = TCcomposite.select('TCG');
  var nbr = TCcomposite.select('NBR');

  // standardize each index - get z-score
  var b5Standard = standardize(b5);
  var tcbStandard = standardize(tcb);
  var tcgStandard = standardize(tcg);
  var nbrStandard = standardize(nbr);

  // combine the standardized TC band collections into a single collection
  var tcStandard = b5Standard.combine(tcbStandard).combine(tcgStandard).combine(nbrStandard);

  // reduce the collection to a single value
  TCcomposite = tcStandard.map(function(img){
    var imgCollection = ee.ImageCollection.fromImages(
      [
        img.select(['B5'],['Z']),//.pow(ee.Image(1)).multiply(img.select('B5').gte(0).where(img.select('B5').lt(0),-1)),
        img.select(['TCB'],['Z']),//.pow(ee.Image(1.5)).multiply(img.select('TCB').gte(0).where(img.select('TCB').lt(0),-1)),
        img.select(['TCG'],['Z']),//.pow(ee.Image(1.5)).multiply(img.select('TCG').gte(0).where(img.select('TCG').lt(0),-1)),
        img.select(['NBR'],['Z'])//.pow(ee.Image(1.5)).multiply(img.select('NBR').gte(0).where(img.select('NBR').lt(0),-1))
      ]
    );
    var reducedImg;
    switch(reducer){
      case 'mean':
        reducedImg = imgCollection.mean();
        break;
      case 'max':
        reducedImg = imgCollection.max();
        break;
      case 'sum':
        reducedImg = imgCollection.sum();
        break;
      default:
        print('The reducer you provided is not supported');
    }

    return reducedImg.multiply(1000).set('system:time_start', img.get('system:time_start'));
  });

  return TCcomposite;
};



// STANDARDIZE A INDEX INDEX - all disturbances are up
var standardizeIndex = function(collection, index){
  var zCollection = collection.map(function(img){
    return calcIndex(img, index, 1);
  });

  zCollection = standardize(zCollection);

  zCollection = zCollection.map(function(img){
    return img.multiply(1000).set('system:time_start', img.get('system:time_start'));
  });

  return zCollection;
};


// BUILD AN LT COLLECTION
var buildLTcollection = function(collection, index, ftvList){
  //print(ftvList)
  var LTcollection;
  switch(index){
    // tasseled cap composite
    case 'TCC':
      LTcollection = makeTCcomposite(collection, 'mean');
      break;
    case 'TCM':
      LTcollection = makeTCcomposite(collection, 'max');
      break;
    case 'TCS':
      LTcollection = makeTCcomposite(collection, 'sum');
      break;

    // 6-band composite - Based on Figure 3 of the linked paper: https://larse.forestry.oregonstate.edu/sites/larse/files/pub_pdfs/Cohen_et_al_2018.pdf
    case 'ENC':
      LTcollection = makeEnsemblecomposite(collection, 'mean');
      break;
    case 'ENM':
      LTcollection = makeEnsemblecomposite(collection, 'max');
      break;
    case 'ENS':
      LTcollection = makeEnsemblecomposite(collection, 'sum');
      break;

    // 6-band composite - Based on Table 5 of the linked paper: https://larse.forestry.oregonstate.edu/sites/larse/files/pub_pdfs/Cohen_et_al_2018.pdf
    case 'ENC1':
      LTcollection = makeEnsemblecomposite1(collection, 'mean');
      break;
    case 'ENM1':
      LTcollection = makeEnsemblecomposite1(collection, 'max');
      break;
    case 'ENS1':
      LTcollection = makeEnsemblecomposite1(collection, 'sum');
      break;

    // standardized versions of indices: mean 0 stdDev 1
    case 'B5z':
      LTcollection = standardizeIndex(collection, 'B5');
      break;
    case 'B7z':
      LTcollection = standardizeIndex(collection, 'B7');
      break;
    case 'TCWz':
      LTcollection = standardizeIndex(collection, 'TCW');
      break;
    case 'TCAz':
      LTcollection = standardizeIndex(collection, 'TCA');
      break;
    case 'NDMIz':
      LTcollection = standardizeIndex(collection, 'NDMI');
      break;
    case 'NBRz':
      LTcollection = standardizeIndex(collection, 'NBR');
      break;

    default:
      //print('default')
      //print(ftvList)
      LTcollection = collection.map(function(img){
              var allStack = calcIndex(img, index, 1);
              var ftvimg;
              for(var ftv in ftvList){
                ftvimg = calcIndex(img, ftvList[ftv], 0)
                                     .select([ftvList[ftv]],['ftv_'+ftvList[ftv].toLowerCase()]);

                  allStack = allStack.addBands(ftvimg)
                                     .set('system:time_start', img.get('system:time_start'));
              }
              return allStack;
      });
    }
  //print(LTcollection)
  return LTcollection;
};

exports.buildLTcollection = buildLTcollection;


// THIS BUILDS AN ANNUAL-SR-TRANSFORMED COLLECTION FOR USE BY TIMESYNC-LEGACY
var buildTSindexCollection = function(collection, ftvList){
  return collection.map(function(img){
            var allStack = ee.Image();
            var ftvimg;
            for(var ftv in ftvList){
              ftvimg = calcIndex(img, ftvList[ftv], 0)
                                   .select([ftvList[ftv]],['ftv_'+ftvList[ftv].toLowerCase()])
                                   .unmask(0);

                allStack = allStack.addBands(ftvimg)
                                   .set('system:time_start', img.get('system:time_start'));
            }
            return allStack.slice(1,null);
         });
};



// TRANSFORM AN ANNUAL SR COLLECTION TO AN ANNUAL COLLECTION OF SELECTED INDICES OR BANDS
var transformSRcollection = function(srCollection, bandList){
  return srCollection.map(function(img){
    var allStack = calcIndex(img, bandList[0], 0);
    for(var band=1; band < bandList.length; band++){
      var bandImg = calcIndex(img, bandList[band], 0);
      allStack = allStack.addBands(bandImg)
                         .set('system:time_start', img.get('system:time_start'));
    }
    return allStack;
  });
};

exports.transformSRcollection = transformSRcollection;



exports.runLT = function(startYear, endYear, startDay, endDay, aoi, index, ftvList, runParams, maskThese){
  maskThese = (typeof maskThese !== 'undefined') ?  maskThese : ['cloud','shadow','snow'];
  var annualSRcollection = buildSRcollection(startYear, endYear, startDay, endDay, aoi, maskThese);
  var annualLTcollection = buildLTcollection(annualSRcollection, index, ftvList);
  runParams.timeSeries = annualLTcollection;
  return ee.Algorithms.TemporalSegmentation.LandTrendr(runParams);
};


// TODO: update the ftv parameter in the guide
var getFittedData = function(lt, startYear, endYear, index, ftv){
  var bandNames = getYearBandNames(startYear, endYear);
  var search;
  if(ftv === true){
    search = index.toUpperCase()+'_fit';
  } else {
    search = 'ftv_'+index.toLowerCase()+'_fit';
  }
  return lt.select(search)
           .arrayFlatten([bandNames]);
};

exports.getFittedData = getFittedData;



// MAKE RGB COMPOSITE
var makeRGBcomposite = function(index, startYear, endYear, startDay, endDay, redYear, greenYear, blueYear, aoi, runParams, nStdev, maskThese){
  var annualSRcollection = buildSRcollection(startYear, endYear, startDay, endDay, aoi, maskThese);
  var annualLTcollection = buildLTcollection(annualSRcollection, index, [index]);
  runParams.timeSeries = annualLTcollection;
  var lt = ee.Algorithms.TemporalSegmentation.LandTrendr(runParams);
  var ftvStack = getFittedData(lt, startYear, endYear, index);
  return ftvStack.select([redYear.toString(),greenYear.toString(),blueYear.toString()]);
};

exports.makeRGBcomposite = makeRGBcomposite;



exports.mapRGBcomposite = function(index, startYear, endYear, startDay, endDay, redYear, greenYear, blueYear, aoi, runParams, nStdev){
  var rgb = makeRGBcomposite(index, startYear, endYear, startDay, endDay, redYear, greenYear, blueYear, aoi, runParams, nStdev);
  var stretch = stdDevStretch(rgb, aoi, nStdev);
  var rgbVis = rgb.visualize({min: stretch[0], max: stretch[1], gamma: 1}).clip(aoi);
  return rgbVis;
  //Map.centerObject(aoi, 10);
  //Map.addLayer(rgbVis);
};


// GENERIC FLATTEN COLLECTION INTO BAND STACK
var collectionToBandStack = function(collection, startYear, endYear, maskFill){
  maskFill = (typeof maskFill !== 'undefined') ?  maskFill : 0;
  var unmaskedCollection = collection.map(function(img){
    return img.unmask(maskFill);
  });
  var collectionArray = unmaskedCollection.toArrayPerBand();
  var nBands = unmaskedCollection.first().bandNames().length().getInfo();
  var years = [];
  var allStack = ee.Image();
  for (var i = startYear; i <= endYear; ++i) years.push(i.toString());
  for (i = 0; i < nBands; i++){
    var bandTS = collectionArray.select([i]).arrayFlatten([years]);
    allStack = ee.Image.cat([allStack, bandTS]);
  }
  return allStack.slice(1,null).toShort();
};
exports.collectionToBandStack = collectionToBandStack;



// FLATTEN SOURCE COLLECTION INTO A GIANT STACK FOR EXPORT - USED BY TIMESYNC
exports.timesyncLegacyStack = function(startYear, endYear, startDay, endDay, aoi, maskThese){
  var annualSRcollection = buildSRcollection(startYear, endYear, startDay, endDay, aoi, maskThese);
  var annualTScollection = buildTSindexCollection(annualSRcollection, ['B1','B2','B3','B4','B5','B7','TCB','TCG','TCW']);

  // THE FOLLOWING COULD BE A GENERIC FUNCTION
  var annualTSArray = annualTScollection.toArrayPerBand();
  var years = [];
  var allStack = ee.Image();
  for (var i = startYear; i <= endYear; ++i) years.push(i.toString());
  for (i = 0; i <= 8; i++){
    var bandTS = annualTSArray.select([i]).arrayFlatten([years]);
    allStack = ee.Image.cat([allStack, bandTS]);
  }

  return allStack.slice(1,null).toShort();
  //////////////////////////////////////////////
};


// #######################################################################################
// ###### PIXEL TIME SERIES FUNCTIONS ####################################################
// #######################################################################################

// RETURN IMAGE BAND VALUES FOR A SINGLE PIXEL AS AN OBJECT
var getPixelInfo = function(img, pixel) {
  return img.unmask(-9999).reduceRegion({
   reducer: 'first',
   geometry: pixel,
   scale: 30
  }).getInfo();
};
exports.getPixelInfo = getPixelInfo;

// RETURN LT RESULTS FOR A SINGLE PIXEL AS AN OBJECT
var ltPixelTimeSeries = function(img, pixel) {
  return img.reduceRegion({
   reducer: 'first',
   geometry: pixel,
   scale: 30
  }).getInfo();
};
exports.ltPixelTimeSeries = ltPixelTimeSeries;

// PARSE OBJECT RETURNED FROM 'getPoint' TO ARRAY OF SOURCE AND FITTED
exports.ltPixelTimeSeriesArray = function(lt, pixel, indexFlip){
  var pixelTS = ltPixelTimeSeries(lt, pixel);
  if(pixelTS.LandTrendr === null){pixelTS.LandTrendr = [[0,0],[0,0],[0,0]]}
  var data = [['Year', 'Original', 'Fitted']];
  var len = pixelTS.LandTrendr[0].length;
  for (var i = 0; i < len; i++) {
    data = data.concat([[pixelTS.LandTrendr[0][i], pixelTS.LandTrendr[1][i]*indexFlip, pixelTS.LandTrendr[2][i]*indexFlip]]);
  }
  return {ts:data, rmse:pixelTS.rmse};
};




// #######################################################################################
// ###### LT ARRAY FUNCTIONS ###############################################################
// #######################################################################################

// get segment info array
var getSegmentData = function(lt, index, delta, right){
  right = (typeof right !== 'undefined') ?  right : false;
  lt = lt.select('LandTrendr');                // select the LandTrendr band
  var vertexMask = lt.arraySlice(0, 3, 4);     // slice out the 'Is Vertex' row - yes(1)/no(0)
  var vertices = lt.arrayMask(vertexMask);     // use the 'Is Vertex' row as a mask for all rows
  var leftList = vertices.arraySlice(1, 0, -1);    // slice out the vertices as the start of segments
  var rightList = vertices.arraySlice(1, 1, null); // slice out the vertices as the end of segments
  var startYear = leftList.arraySlice(0, 0, 1);    // get year dimension of LT data from the segment start vertices
  var startVal = leftList.arraySlice(0, 2, 3);     // get spectral index dimension of LT data from the segment start vertices
  var endYear = rightList.arraySlice(0, 0, 1);     // get year dimension of LT data from the segment end vertices
  var endVal = rightList.arraySlice(0, 2, 3);      // get spectral index dimension of LT data from the segment end vertices
  var dur = endYear.subtract(startYear);       // subtract the segment start year from the segment end year to calculate the duration of segments
  var mag = endVal.subtract(startVal);         // substract the segment start index value from the segment end index value to calculate the delta of segments
  var rate = mag.divide(dur);                  // calculate the rate of spectral change

  if(delta.toLowerCase() == 'all'){
    if((indexFlipper(index) == -1) && (right === true)){
      startVal = startVal.multiply(-1);
      endVal = endVal.multiply(-1);
      mag = mag.multiply(-1);
      rate = rate.multiply(-1);
    }
    return ee.Image.cat([startYear.add(1), endYear, startVal, endVal, mag, dur, rate])
                   .toArray(0);
  } else if((delta.toLowerCase() == 'gain') || (delta.toLowerCase() == 'loss')){
    var changeTypeMask;
    if(delta.toLowerCase() == 'gain'){
       changeTypeMask = mag.lt(0);
    } else if(delta.toLowerCase() == 'loss'){
      changeTypeMask = mag.gt(0);
    }

    var flip = indexFlipper(index);
    return ee.Image.cat([
      startYear.add(1).arrayMask(changeTypeMask).unmask(ee.Image(ee.Array([[-9999]]))),
      endYear.arrayMask(changeTypeMask).unmask(ee.Image(ee.Array([[-9999]]))),
      startVal.arrayMask(changeTypeMask).multiply(flip).unmask(ee.Image(ee.Array([[-9999]]))),
      endVal.arrayMask(changeTypeMask).multiply(flip).unmask(ee.Image(ee.Array([[-9999]]))),
      mag.arrayMask(changeTypeMask).abs().unmask(ee.Image(ee.Array([[-9999]]))),
      dur.arrayMask(changeTypeMask).unmask(ee.Image(ee.Array([[-9999]]))),
      rate.arrayMask(changeTypeMask).abs().unmask(ee.Image(ee.Array([[-9999]])))
    ]).toArray(0);
  }
};
exports.getSegmentData = getSegmentData;



// #######################################################################################
// ###### HELPER FUNCTIONS ###############################################################
// #######################################################################################

// GET A SERIES OF BANDS NAMES AS YEARS
var getYearBandNames = function(startYear, endYear){
  var years = [];                                                           // make an empty array to hold year band names
  for (var i = startYear; i <= endYear; ++i) years.push(i.toString()); // fill the array with years from the startYear to the endYear and convert them to string
  return years;
};
exports.getYearBandNames = getYearBandNames;


// STANDARD DEVIATION STRETCH
var stdDevStretch = function(img, aoi, nStdev){
  var mean = img.reduceRegion({reducer:ee.Reducer.mean(), geometry:aoi, scale:900, bestEffort:true, maxPixels:1e4})
              .toArray()
              .reduce(ee.Reducer.mean(), [0]);

  var stdDev = img.reduceRegion({reducer:ee.Reducer.stdDev(), geometry:aoi, scale:900, bestEffort:true, maxPixels:1e4})
                .toArray()
                .reduce(ee.Reducer.mean(), [0])
                .multiply(nStdev);

  var max = mean.add(stdDev).getInfo()[0];
  var min = mean.subtract(stdDev).getInfo()[0];
  return [min, max];
};

// INDEX FLIPPER
var indexFlipper = function(index){
  var indexObj = {'NBR':-1,'NDVI':-1,'EVI':-1, 'NDSI':-1,
                  'TCB':1,'TCG':-1,'TCW':-1,
                  'B1':1,'B2':1,'B3':1,'B4':-1,'B5':1,'B7':1,
                  'ENC':1,'ENC1':1,'TCC':1,
                  'NBRz':1,'B5z':1};
  return indexObj[index];
};
exports.indexFlipper = indexFlipper;


var makeBoolean = function(value){
  if(typeof(value) !== "boolean"){
    value = value.toLowerCase() != 'false';
  }
  return value;
};
exports.makeBoolean = makeBoolean;
