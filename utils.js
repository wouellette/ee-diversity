// ****************************************************************************************************************** //
// *********************************** Utilities for the ee-diversity package *************************************** //
// ****************************************************************************************************************** //

var PI = ee.Number(3.14159265359);
var MAX_SATELLITE_ZENITH = 7.5;
var MAX_DISTANCE = 1000000;
var UPPER_LEFT = 0;
var LOWER_LEFT = 1;
var LOWER_RIGHT = 2;
var UPPER_RIGHT = 3;

// Function to load Sentinel-2 data and its corresponding cloud probability information, based on an area and time range
exports.loadImageCollection = function(collection_name, time_range, cloud_filter, geom){
  // Import Sentinel-2 L1C or L2A, depending on year selected.
  var s2 = ee.ImageCollection(collection_name)
    .filterDate(time_range.get('start'), time_range.get('end'))
    .filterBounds(geom)
    .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', cloud_filter);

  // Import and filter s2cloudless.
  var s2_cloudless_col = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
      .filterBounds(geom)
      .filterDate(time_range.get('start'), time_range.get('end'));

  // Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
  var s2_cl = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
      primary: s2,
      secondary: s2_cloudless_col,
      condition: ee.Filter.equals({
          leftField: 'system:index',
          rightField: 'system:index'})
  })).sort('system:time_start');

  return s2_cl
}

/**
 * Generate a Sentinel-2 Cloud Mask Time Series using the Sentinelhub Cloud Mask
 * @param {ImageCollection} collection: The input Sentinel-2 Collection
 * @param {List} bands: The bands to select from the Sentinel-2 collection
 * @param {String} start_year: start year of the Sentinel-2 time series
 * @param {String} end_year: end year of the Sentinel-2 time series
 * @param {FeatureCollection} geom: Input geometry to stratify the time series
 * @returns {ImageCollection}: A cloud-masked Sentinel-2 time series
 * @ignore
*/
exports.S2CloudMasked = function(collection, bands, start_year, end_year, surface_water,
                                 cloud_filter, nir_drk_thresh, cld_prj_dist, cld_prb_thresh, cld_buffer,
                                 snow_thresh, mask_res, geom){
    var date_range = ee.Dictionary({'start': start_year+'-01-01', 'end': end_year+'-12-31'});

    // Load sentinel-2 data for the full year of interest
    var s2_cl = exports.loadImageCollection(collection, date_range, cloud_filter, geom.geometry());

    // Perform cloud masking using the S2 cloud probabilities assets from s2cloudless,
    // courtesy of Sentinelhub/EU/Copernicus/ESA
    var masked_collection = s2_cl.filter(ee.Filter.notNull(['MEAN_INCIDENCE_AZIMUTH_ANGLE_B3',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B4',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B5',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B6',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B7',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B8A',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B11',
                                                            'MEAN_INCIDENCE_AZIMUTH_ANGLE_B12',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B3',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B4',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B5',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B6',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B7',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B8A',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B11',
                                                            'MEAN_INCIDENCE_ZENITH_ANGLE_B12',
                                                            'MEAN_SOLAR_AZIMUTH_ANGLE',
                                                            'MEAN_SOLAR_ZENITH_ANGLE']))
      .map(exports.addCloudShadowMask(surface_water.not(), nir_drk_thresh, cld_prj_dist,
                                      cld_prb_thresh, cld_buffer, snow_thresh, mask_res, 1e4))
      .map(exports.applyCloudShadowMask)
      .select(bands);

    return masked_collection
};

// Function to combine shadow and cloud masks to image.
exports.addCloudShadowMask = function(water_valmask,  nir_drk_thresh, cld_prj_dist,
                                      cld_prb_thresh, cld_buffer, snow_thresh, mask_res, sr_band_scale){
    // water_valmask : water validity mask, indicating locations of non-water pixels for the cloud shadow detection
    // sr_band_scale: scaling factor. 10000 for Sentinel-2 GEE assets

    sr_band_scale = sr_band_scale || 1;

    var wrap = function(img){
      // img: A sentinel-2 image

      // Add cloud component bands.
      var img_cloud = _addCloudBands(img, cld_prb_thresh);

      // Add cloud shadow component bands.
      var img_cloud_shadow = _addShadowBands(img_cloud, water_valmask,
                                             nir_drk_thresh, cld_prj_dist, mask_res, sr_band_scale);

      var img_snow = _addSnowBands(img_cloud_shadow, snow_thresh);

      // Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
      var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0);

      // Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
      // 60 m scale is for speed, and assumes clouds don't require 10 m precision.
      is_cld_shdw = is_cld_shdw.focal_min(2) // Morphological Opening operation is an erosion (focal_min) followed by
                    .focal_max(cld_buffer * 2 / mask_res) // a dilation (focal_max)
                    .reproject({crs: img.select([0]).projection(), scale: mask_res}) // reproject to resample resolution
                    .rename('cloudmask');

      var is_snw = img_snow.select('snow').rename('snowmask');

      // Add the final cloud-shadow mask to the image.
      return img_cloud_shadow.addBands(is_cld_shdw).addBands(is_snw)
      }
    return wrap
  };


 // Function to apply the final cloud mask to the image.
exports.applyCloudShadowMask = function(img){
  // img: A sentinel-2 image

  // Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
  var not_cld_shdw = img.select('cloudmask').not();
  var not_snw = img.select('snowmask').not();

  // Subset reflectance bands and update their masks, return the result.
  return img.updateMask(not_cld_shdw.and(not_snw))
  }

// Function to add the cloud probability band to the image based on the specified cloud probability threshold.
function _addCloudBands(img, cld_prb_thresh){
    // img: A sentinel-2 image

    // Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img.get('s2cloudless')).select('probability');

    // Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(cld_prb_thresh).rename('clouds');

    // Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))
}

// Function to add the cloud shadow mask to the image based on specified parameters.
function _addShadowBands(img, water_valmask, nir_drk_thresh, cld_prj_dist, mask_res, sr_band_scale){
    // img: A sentinel-2 image
    // water_valmask : water validity mask, indicating locations of non-water pixels for the cloud shadow detection
    // sr_band_scale: scaling factor. 10000 for Sentinel-2 GEE assets

    // Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var dark_pixels = img.select('B8').lt(nir_drk_thresh*sr_band_scale)
                      .multiply(water_valmask).rename('dark_pixels');

    // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = img.select('clouds').directionalDistanceTransform(shadow_azimuth, cld_prj_dist / 10) // Looks in az
                                                                                                        // direction
                   .reproject({crs:img.select(0).projection(), scale: mask_res}) // reproject to reduce proc time
                   .select('distance').mask().rename('cloud_transform');

    // Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows');

    // Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
}

function _addSnowBands(img, snow_thresh){
  return img.addBands(img.normalizedDifference(['B3', 'B11']).gt(snow_thresh).rename('snow'));
}

exports.applyBRDF = function(col_flag){

  var wrap = function(image){
    var date = image.date();
    var footprint = ee.List(image.geometry().bounds().coordinates().get(0));
    var angles =  getsunAngles(date, footprint);
    var sunAz = angles[0];
    var sunZen = angles[1];

    var viewAz = azimuth(footprint);
    var viewZen = zenith(footprint);

    var kval = _kvol(sunAz, sunZen, viewAz, viewZen);
    var kvol = kval[0];
    var kvol0 = kval[1];
    var result = _apply(image, kvol.multiply(PI), kvol0.multiply(PI), col_flag);

    return result;
  }

  return wrap;
}

/* Get sunAngles from the map given the data.
*
* date:  ee.date object
* footprint: geometry of the image
*/
function getsunAngles(date, footprint){
  var jdp = date.getFraction('year');
  var seconds_in_hour = 3600;
  var  hourGMT = ee.Number(date.getRelative('second', 'day')).divide(seconds_in_hour);

  var latRad = ee.Image.pixelLonLat().select('latitude').multiply(PI.divide(180));
  var longDeg = ee.Image.pixelLonLat().select('longitude');

  // Julian day proportion in radians
  var jdpr = jdp.multiply(PI).multiply(2);

  var a = ee.List([0.000075, 0.001868, 0.032077, 0.014615, 0.040849]);
  var meanSolarTime = longDeg.divide(15.0).add(ee.Number(hourGMT));
  var localSolarDiff1 = value(a, 0)
          .add(value(a, 1).multiply(jdpr.cos()))
          .subtract(value(a, 2).multiply(jdpr.sin()))
          .subtract(value(a, 3).multiply(jdpr.multiply(2).cos()))
          .subtract(value(a, 4).multiply(jdpr.multiply(2).sin()));

  var localSolarDiff2 = localSolarDiff1.multiply(12 * 60);

  var localSolarDiff = localSolarDiff2.divide(PI);
  var trueSolarTime = meanSolarTime
          .add(localSolarDiff.divide(60))
          .subtract(12.0);

  // Hour as an angle;
  var ah = trueSolarTime.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2).multiply(PI.divide(180))) ;
  var b = ee.List([0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480]);
  var delta = value(b, 0)
        .subtract(value(b, 1).multiply(jdpr.cos()))
        .add(value(b, 2).multiply(jdpr.sin()))
        .subtract(value(b, 3).multiply(jdpr.multiply(2).cos()))
        .add(value(b, 4).multiply(jdpr.multiply(2).sin()))
        .subtract(value(b, 5).multiply(jdpr.multiply(3).cos()))
        .add(value(b, 6).multiply(jdpr.multiply(3).sin()));

  var cosSunZen = latRad.sin().multiply(delta.sin())
        .add(latRad.cos().multiply(ah.cos()).multiply(delta.cos()));
  var sunZen = cosSunZen.acos();

  // sun azimuth from south, turning west
  var sinSunAzSW = ah.sin().multiply(delta.cos()).divide(sunZen.sin());
  sinSunAzSW = sinSunAzSW.clamp(-1.0, 1.0);

  var cosSunAzSW = (latRad.cos().multiply(-1).multiply(delta.sin())
                    .add(latRad.sin().multiply(delta.cos()).multiply(ah.cos())))
                    .divide(sunZen.sin());
  var sunAzSW = sinSunAzSW.asin();

  sunAzSW = where(cosSunAzSW.lte(0), sunAzSW.multiply(-1).add(PI), sunAzSW);
  sunAzSW = where(cosSunAzSW.gt(0).and(sinSunAzSW.lte(0)), sunAzSW.add(PI.multiply(2)), sunAzSW);

  var sunAz = sunAzSW.add(PI);
   // # Keep within [0, 2pi] range
    sunAz = where(sunAz.gt(PI.multiply(2)), sunAz.subtract(PI.multiply(2)), sunAz);

  var footprint_polygon = ee.Geometry.Polygon(footprint);
  sunAz = sunAz.clip(footprint_polygon);
  sunAz = sunAz.rename(['sunAz']);
  sunZen = sunZen.clip(footprint_polygon).rename(['sunZen']);

  return [sunAz, sunZen];
}


/* Get azimuth.
*
*
* footprint: geometry of the image
*/
function azimuth(footprint){
    function x(point){return ee.Number(ee.List(point).get(0))}
    function  y(point){return ee.Number(ee.List(point).get(1))}

    var upperCenter = line_from_coords(footprint, UPPER_LEFT, UPPER_RIGHT).centroid(0.001).coordinates();
    var lowerCenter = line_from_coords(footprint, LOWER_LEFT, LOWER_RIGHT).centroid(0.001).coordinates();
    var slope = ((y(lowerCenter)).subtract(y(upperCenter))).divide((x(lowerCenter)).subtract(x(upperCenter)));
    var slopePerp = ee.Number(-1).divide(slope);
    var azimuthLeft = ee.Image(PI.divide(2).subtract((slopePerp).atan()));
    return azimuthLeft.rename(['viewAz']);
  }

/* Get zenith.
*
*
* footprint: geometry of the image
*/
function zenith(footprint){
    var leftLine = line_from_coords(footprint, UPPER_LEFT, LOWER_LEFT);
    var rightLine = line_from_coords(footprint, UPPER_RIGHT, LOWER_RIGHT);
    var leftDistance = ee.FeatureCollection(leftLine).distance(MAX_DISTANCE);
    var rightDistance = ee.FeatureCollection(rightLine).distance(MAX_DISTANCE);
    var viewZenith = rightDistance.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2))
          .divide(rightDistance.add(leftDistance))
          .subtract(ee.Number(MAX_SATELLITE_ZENITH))
          .clip(ee.Geometry.Polygon(footprint))
          .rename(['viewZen']);
    return viewZenith.multiply(PI.divide(180));
}

  /* apply function to all bands
  *
  * http://www.mdpi.com/2072-4292/9/12/1325/htm#sec3dot2-remotesensing-09-01325
  * https://www.sciencedirect.com/science/article/pii/S0034425717302791
  *
  * image : the image to apply the function to
  * kvol:
  * kvol0
  *
  */
function _apply(image, kvol, kvol0, col_flag){
      var f_iso = 0;
      var f_geo = 0;
      var f_vol = 0;
			var blue = _correct_band(image, 'B2', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372);
			var green = _correct_band(image, 'B3', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580);
			var red = _correct_band(image, 'B4', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574);
			var re1 = ee.Image(ee.Algorithms.If(ee.String(col_flag).equals('s2'), _correct_band(image, 'B5', kvol, kvol0, f_iso=0.2085, f_geo=0.0256, f_vol=0.0845), ee.Image(0)));
			var re2 = ee.Image(ee.Algorithms.If(ee.String(col_flag).equals('s2'), _correct_band(image, 'B6', kvol, kvol0, f_iso=0.2316, f_geo=0.0273, f_vol=0.1003), ee.Image(0)));
			var re3 = ee.Image(ee.Algorithms.If(ee.String(col_flag).equals('s2'), _correct_band(image, 'B7', kvol, kvol0, f_iso=0.2599, f_geo=0.0294, f_vol=0.1197), ee.Image(0)));
      var nir = _correct_band(image, 'B8', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535);
      var re4 = ee.Image(ee.Algorithms.If(ee.String(col_flag).equals('s2'), _correct_band(image, 'B8A', kvol, kvol0, f_iso=0.2907, f_geo=0.0410, f_vol=0.1611), ee.Image(0)));
      var swir1 = _correct_band(image, 'B11', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154);
      var swir2 = _correct_band(image, 'B12', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639);
      var brdf = ee.Image(ee.Algorithms.If(ee.String(col_flag).equals('s2'),
                          image.addBands({srcImg: ee.Image([blue, green, red, re1, re2, re3, nir, re4, swir1, swir2]), overwrite: true}),
                          image.addBands({srcImg: ee.Image([blue, green, red, nir, swir1, swir2]), overwrite: true})));
			return brdf
}

/* correct band function
  *
  *
  * image : the image to apply the function to
  * band_name
  * kvol
  * kvol0
  * f_iso
  * f_geo
  * f_vol
  *
  */
function _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol){
	//"""fiso + fvol * kvol + fgeo * kgeo"""
	var iso = ee.Image(f_iso);
	var geo = ee.Image(f_geo);
	var vol = ee.Image(f_vol);
	var pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred']);
	var pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0']);
	var cfac = pred0.divide(pred).rename(['cfac']);
	var corr = image.select(band_name).multiply(cfac).rename([band_name]);
	return corr;
}

/* calculate kvol and kvol0
*
* sunAZ
* sunZen
* viewAz
* viewZen
*
*/
function _kvol(sunAz, sunZen, viewAz, viewZen){
	//"""Calculate kvol kernel.
	//From Lucht et al. 2000
	//Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""

	var relative_azimuth = sunAz.subtract(viewAz).rename(['relAz']);
	var pa1 = viewZen.cos().multiply(sunZen.cos());
	var pa2 = viewZen.sin().multiply(sunZen.sin()).multiply(relative_azimuth.cos());
	var phase_angle1 = pa1.add(pa2);
	var phase_angle = phase_angle1.acos();
	var p1 = ee.Image(PI.divide(2)).subtract(phase_angle);
	var p2 = p1.multiply(phase_angle1);
	var p3 = p2.add(phase_angle.sin());
	var p4 = sunZen.cos().add(viewZen.cos());
	var p5 = ee.Image(PI.divide(4));

	var kvol = p3.divide(p4).subtract(p5).rename(['kvol']);

	var viewZen0 = ee.Image(0);
	var pa10 = viewZen0.cos().multiply(sunZen.cos());
	var pa20 = viewZen0.sin().multiply(sunZen.sin()).multiply(relative_azimuth.cos());
	var phase_angle10 = pa10.add(pa20);
	var phase_angle0 = phase_angle10.acos();
	var p10 = ee.Image(PI.divide(2)).subtract(phase_angle0);
	var p20 = p10.multiply(phase_angle10);
	var p30 = p20.add(phase_angle0.sin());
	var p40 = sunZen.cos().add(viewZen0.cos());
	var p50 = ee.Image(PI.divide(4));

	var kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0']);

	return [kvol, kvol0]}

/* helper function
*
*
*
*/
function line_from_coords(coordinates, fromIndex, toIndex){
    return ee.Geometry.LineString(ee.List([
      coordinates.get(fromIndex),
      coordinates.get(toIndex)]));
}

function where(condition, trueValue, falseValue){
  var trueMasked = trueValue.mask(condition);
  var falseMasked = falseValue.mask(invertMask(condition));
      return trueMasked.unmask(falseMasked);
}

function invertMask(mask){
    return mask.multiply(-1).add(1);
}


function value(list,index){
    return ee.Number(list.get(index));
}

/**
 * Function to compute trend (i.e. slope of change) from a time series.
   This can be used to compute the productivity trend from SDG 15.3.1, but can be used to derive others trends too
 * @param {ImageCollection} collection: The input time series used to compute the SeRGS
 * @param {String} independent_var: The independent variable to use
 * @param {String} dependent_var: The dependent variable to use
 * @returns {Image}: The trend output
 * @ignore
*/
exports.trendTS = function(collection, independent_var, dependent_var){

    var trend_mk = exports.mannKendall(collection.select(dependent_var)).rename('mannkendall');

    var trend_ols = collection.select([independent_var, dependent_var])
                              .reduce(ee.Reducer.linearFit(), 16)
                              .select('scale')
                              .addBands(trend_mk);

    var trend_signif = exports.signifMask(trend_ols.select('scale'), trend_mk, collection.size())
                            .rename('significance');

    trend_ols = trend_ols.addBands(trend_signif);

    return trend_ols;
}

/**
 * Calculate Mann Kendall's S statistic.
 This function returns the Mann Kendall's S statistic, assuming that n is
 less than 40. The significance of a calculated S statistic is found in
 table A.30 of Nonparametric Statistical Methods, second edition by
 Hollander & Wolfe.
 * @param {imageCollection} collection: Input image collection for which to calculate the mann-kendall S statistic.
 * @returns {Image} mk_stat: The Mann Kendall S statistic.
*/
exports.mannKendall = function(collection){
  var afterFilter = ee.Filter.lessThan({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

  var joined = ee.ImageCollection(ee.Join.saveAll('after').apply({
    primary: collection,
    secondary: collection,
    condition: afterFilter
  }));

  var sign = function(i, j) { // i and j are images
    //return ee.Image(j).neq(i) // Zero case
    //    .multiply(ee.Image(j).subtract(i).clamp(-1, 1)).int();
      var concordant = ee.Image(i).lt(j).rename('concordant');
      var discordant = ee.Image(i).gt(j).rename('discordant');
      return concordant.addBands(discordant);
  };

  var mk = ee.ImageCollection(joined.map(function(current) {
    var afterCollection = ee.ImageCollection.fromImages(current.get('after'));
    return afterCollection.map(function(image) {
      // The unmask is to prevent accumulation of masked pixels that
      // result from the undefined case of when either current or image
      // is masked.  It won't affect the sum, since it's unmasked to zero.
      return ee.Image(sign(current, image)).unmask(0);
    });
    // Set parallelScale to avoid User memory limit exceeded.
  }).flatten()).reduce('sum', 4);

  var mk_stat = mk.select('concordant_sum').subtract(mk.select('discordant_sum'));
  return mk_stat.toFloat()
}

/**
 * Generate a significance mask from the mann-kendall test results
 * @param {Image} img: Input image for which to generate a statistical significance mask.
 * @param {Image} mk_trend: Mann-Kendall S statistic image.
 * @returns {Image}: The Statistical significance mask for the 90%, 95% and 99% confidence intervals.
*/
exports.signifMask = function(img, mk_trend, period){
  // Define Kendall parameter values for a significance of 0.05
  //var period = end_year.get('year').subtract(start_year.get('year')).add(1);
  var kendall90 = ee.Number(_kendallCoefficient(period, 90));
  var kendall95 = ee.Number(_kendallCoefficient(period, 95));
  var kendall99 = ee.Number(_kendallCoefficient(period, 99));
  // Create final productivity trajectory output layer. Positive values are
  // significant increase, negative values are significant decrease.
  return ee.Image(-32768)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall90)), 1)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall95)), 2)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall99)), 3)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall90)), -1)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall95)), -2)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall99)), -3)
        .where(mk_trend.abs().lte(kendall90), 0);
}

/**
 * Hard-coded Kendall Coefficients as look-up dictionary
 * @param {Image} n: The kendall coefficient number to look up.
 * @param {Image} level: The confidence interval to use (90%, 95% or 99%).
 * @returns {Image}: The looked-up coefficient value.
*/
function _kendallCoefficient(n, level){
    // The minus 4 is because the indexing below for a sample size of 4
    n = n.subtract(4);
    var coefs = {90: ee.List([4, 6, 7, 9, 10, 12, 15, 17, 18, 22, 23, 27, 28, 32, 35, 37, 40, 42,
                  45, 49, 52, 56, 59, 61, 66, 68, 73, 75, 80, 84, 87, 91, 94, 98, 103,
                  107, 110, 114, 119, 123, 128, 132, 135, 141, 144, 150, 153, 159,
                  162, 168, 173, 177, 182, 186, 191, 197, 202]),
               95: ee.List([4, 6, 9, 11, 14, 16, 19, 21, 24, 26, 31, 33, 36, 40, 43, 47, 50, 54,
                    59, 63, 66, 70, 75, 79, 84, 88, 93, 97, 102, 106, 111, 115, 120,
                    126, 131, 137, 142, 146, 151, 157, 162, 168, 173, 179, 186, 190,
                    197, 203, 208, 214, 221, 227, 232, 240, 245, 251, 258]),
               99: ee.List([6, 8, 11, 18, 22, 25, 29, 34, 38, 41, 47, 50, 56, 61, 65, 70, 76, 81,
                    87, 92, 98, 105, 111, 116, 124, 129, 135, 142, 150, 155, 163, 170,
                    176, 183, 191, 198, 206, 213, 221, 228, 236, 245, 253, 260, 268,
                    277, 285, 294, 302, 311, 319, 328, 336, 345, 355, 364])}
    return coefs[level].get(n);
}

/**
 * Load Annual Land Cover data from WorldCover, Dynamic World and GLAD Land Cover Change 2000-2020
 * @param {String} year: The year for which to load the land cover data
 * @param {FeatureCollection} geom: Geometry to stratify the land cover datasets
 * @returns {Array} landcover_data: an array containing World Cover-related data,
                                    Dynamic World-related data, and GLAD Land Cover Change data.
 * @ignore
 */
exports.landCoverDatasets = function(year, geom){

    // World Cover Dataset
    var world_cover = _worldCover(year);

    // Dynamic World Cover Dataset
    var dw = _dynamicWorld(year, geom);
    var dynamic_world = dw[0];
    var proba_hillshade = dw[1];

    // GLAD land cover change dataset 2000-2020
    var landmask = ee.Image("projects/glad/landBuffer4").mask();
    var LCLUC2020 = ee.Image('projects/glad/GLCLU2020/LCLUC_2020').updateMask(landmask);
    var wetland_cover = world_cover.eq(7).or(LCLUC2020.gte(100).and(LCLUC2020.lt(150))).selfMask();
    var tree_cover = LCLUC2020.gte(25).and(LCLUC2020.lte(100))
                     .or(LCLUC2020.gte(125).and(LCLUC2020.lte(200)))
                     .selfMask().rename('constant');

    var glad_lc_change = _gladLCChange(landmask);

    // Conflate the World Cover 2021 data with GLAD LULC2020 permanent water extent, which is more complete.
    world_cover = world_cover.where(LCLUC2020.gt(100).and(LCLUC2020.lt(150)), 6);

    var forest_change = ee.Image("UMD/hansen/global_forest_change_2020_v1_8");

    var landcover_data = [world_cover, tree_cover, wetland_cover,
                          dynamic_world, proba_hillshade, glad_lc_change, forest_change];

    return landcover_data
}

/**
 * Load Dynamic World Cover Trends from Land Cover Class Probabilities
 * @param {String} start_year: The start year for which to load the land cover data
 * @param {String} end_year: The end year for which to load the land cover data
 * @param {String} agg_interval: Aggregation interval in number of days for each trend timestep
 * @param {FeatureCollection} geom: Geometry to stratify the land cover datasets
 * @returns {Dictionary} dw_ols: Dictionary containing the class-specific trends (slope of change)
                                 with associated significance masks.
 * @ignore
 */
exports.dynamicWorldTrend = function(start_year, end_year, agg_interval, geom){
    // Generate temporal land cover trends from Dynamic World probabilities
    var dw_class_names = ['water', 'trees', 'grass', 'flooded_vegetation', 'crops',
                                    'shrub_and_scrub', 'built', 'bare', 'snow_and_ice'];

    var dw_col = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
                  .filterBounds(geom.geometry())
                  .filterDate(start_year+'-01-01', end_year+'-12-31T23:59:59');

    // Generate temporal intervals based on the AGG_INTERVAL provided
    // (default is 365/6 = 6 time intervals of 2 months each for a year)
    var time_intervals = exports.extractTimeRanges(start_year+'-01-01', end_year+'-12-31', agg_interval);
    // Generate harmonized monthly time series of FCover as input to the
    // vegetation factor V of the sustainability factor S of the RUSLE equation
    var dw_ts = ARDS2.harmonizedTS(dw_col,
                                        dw_class_names,
                                        time_intervals,
                                        {agg_type: 'geomedian'})
                .map(exports.createTimeBand);

    // Dynamic World Probability Trends
    var mk_trend_dw;
    var dw_trend;
    var signif_dw;
    var dw_ols = {};
    for (var i = 0; i < dw_class_names.length; i++) {
      dw_trend = landDegradation.trendTS(dw_ts, 'system:time_start', dw_class_names[i]);
      dw_ols[dw_class_names[i]] = dw_trend;
    }

    return dw_ols;
}

/**
 * Nested function to load the World Cover data
 * @param {String} year: The year for which to load the land cover data
 * @returns {Image} world_cover: World Cover output
 * @ignore
*/
function _worldCover(year){
    var world_cover;
    if (year > '2020') {
        world_cover = ee.Image("ESA/WorldCover/v200/2021")
                      //.filterBounds(country.geometry()).mosaic()
                      // Remapping from world cover classes to a harmonized land cover nomenclature
                      .remap([10,20,30,40,50,60,70,80,90,95,100],
                             [3,4,5,2,1,8,9,6,7,3,5])
                      .rename('constant');
    } else {
        world_cover = ee.Image("ESA/WorldCover/v200/2020")
                      //.filterBounds(country.geometry()).mosaic()
                      // Remapping from world cover classes to a harmonized land cover nomenclature
                      .remap([10,20,30,40,50,60,70,80,90,95,100],
                             [3,4,5,2,1,8,9,6,7,3,5])
                      .rename('constant');
    }

    return world_cover
}

/**
 * Nested function to load the Dynamic World Data, aggregated for a specific year
 * @param {String} year: The year for which to load the land cover data
 * @param {FeatureCollection} geom: Geometry to stratify the land cover datasets
 * @returns {Image} world_cover: World Cover output
 * @ignore
*/
function _dynamicWorld(year, geom){
    // Alternative to WorldCover, which provides a yearly land cover dataset
    // (to be further tested whether it is applicable across geographies)
    // Temporal aggregation process tailored for the Sahelian context

    var dw_class_names = ['water', 'trees', 'grass', 'flooded_vegetation', 'crops',
                                    'shrub_and_scrub', 'built', 'bare', 'snow_and_ice'];

    var dw_col = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
                  .filterBounds(geom.geometry())
                  .filterDate(year+'-01-01', year+'-12-31T23:59:59');

    var dw_col_filt = dw_col.map(function(img){
      var img_mask = img.select(dw_class_names).updateMask(img.select(dw_class_names).gte(0.4));
      var img_mask_bs = img.select('bare').updateMask(img.select('bare').gte(0.6));
      img_mask = img_mask.addBands(img_mask_bs, ['bare'], true);
      return img_mask.addBands(img.select('label').updateMask(img_mask.reduce(ee.Reducer.max()).gt(0)));
    });

    // Get the most likely class probability.
    var top1Prob = dw_col.select(dw_class_names).reduce(ee.Reducer.mean(), 4).regexpRename('_mean', '');
    var top1Tree = dw_col_filt.select('trees').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');
    var top1Grass = dw_col_filt.select('grass').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');
    var top1Shrub = dw_col_filt.select('shrub_and_scrub').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');
    var top1Label = dw_col_filt.select('label').reduce(ee.Reducer.mode(), 4).regexpRename('_mode', '');
    var top1All = dw_col.select('label').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');

    // Create a hillshade of the most likely class probability on [0, 1];
    var proba_hillshade = ee.Terrain.hillshade(top1Prob.reduce(ee.Reducer.max()).multiply(100)
                                                 .reproject({crs:ee.Image(dw_col.first()).projection(), scale:10}))
                                                 .rename('proba');

    // Create an RGB image of the label (most likely class) on [0, 1].
    var dynamic_world = top1Label.remap([0, 1, 2, 3, 4, 5, 6, 7, 8], [6, 3, 5, 7, 2, 4, 1, 8, 9]).rename('label')
                       .where(top1Tree.divide(top1All).lt(0.6)
                              .and(top1Prob.select('crops').gt(top1Prob.select('shrub_and_scrub')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('grass')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('bare')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('water')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('flooded_vegetation')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('built')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('snow_and_ice'))), 2)
                       .where(top1Tree.divide(top1All).lt(0.6)
                              .and(top1Prob.select('crops').subtract(top1Prob.select('bare')).abs().lt(0.2))
                              .and(top1Prob.select('shrub_and_scrub').gt(top1Prob.select('grass')))
                              //.and(top1Prob.select('shrub_and_scrub').gt(top1Prob.select('crops'))))
                              .and(top1Shrub.gte(1))
                              .and(top1Grass.lt(top1Shrub)), 4)
                       .where(top1Tree.divide(top1All).lt(0.6)
                              .and(top1Prob.select('crops').subtract(top1Prob.select('bare')).abs().lt(0.2))
                              .and(top1Prob.select('grass').gt(top1Prob.select('shrub_and_scrub')))
                              //.and(top1Prob.select('grass').gt(top1Prob.select('crops')))
                              .and(top1Grass.gte(1))
                              .and(top1Shrub.lt(top1Grass)), 5)
                       .unmask(8);

    return [dynamic_world, proba_hillshade]
}

/**
 * Nested function to load the GLAD Land Cover Change
 * @param {String} landmask: A land surface mask to mask out the coastal water
 * @returns {Image} change_map: Reclassified GLAD Change data output
 * @ignore
*/
function _gladLCChange(landmask){

    // This data is used to refine land cover data, especially wrt the water extents
    var change = ee.Image('projects/glad/GLCLU2020/LCLUC').updateMask(landmask).rename('constant');

    var change_map = change.gt(50).and(change.lt(100)).multiply(7) // forest gain
                   .add(change.gt(150).and(change.lt(200)).multiply(8)) // wetland gain
                   .add(change.eq(240).multiply(1)) // forest --> short vegetation
                   .add(change.eq(247).multiply(4)) // short vegetation --> cropland
                   .add(change.eq(245).multiply(2)) // forest --> cropland
                   .add(change.eq(251).multiply(6)) // built-up gain
                   .add(change.eq(249).multiply(5)) // cropland --> short vegetation
                   .add(change.eq(246).multiply(3)) // wetland --> cropland
                   .add(change.eq(209).multiply(9)) // water loss
                   .add(change.eq(210).multiply(10)) // water gain
                   .add(change.eq(211).multiply(11)) // Ephemeral water
                   .selfMask()
                   .rename('constant');

    return change_map;
}

/**
 * This function adds a time band to the image
 * @param {Image} img: Input image for which to create time band
 * @returns {Image}: The original image conflated with a time band expressed in years
 * @ignore
*/
exports.createTimeBand = function(img) {
  // Scale milliseconds by a large constant to avoid very small slopes
  // in the linear regression output.
  return img.addBands(img.metadata('system:time_start').divide(3.154e10).add(1970));
};

//Define function to extract time intervals to use to generate the temporal composites from Sentinel collections
exports.extractTimeRanges = function(start, end, agg_interval){
    /*
    Extract the time range data from the received time range and aggregation interval e.g.,
    input time interval: time_interval = ['2019-01-01','2020-01-01'], agg_interval: 60 days
    generate the following time intervals:
    time_range = [("2019-01-01T00:00:00Z", "2019-03-01T00:00:00Z"),
                ("2019-03-01T00:00:00Z", "2019-05-01T00:00:00Z"),
                ("2019-05-01T00:00:00Z", "2019-07-01T00:00:00Z"),
                ("2019-07-01T00:00:00Z", "2019-09-01T00:00:00Z"),
                ("2019-09-01T00:00:00Z", "2019-11-01T00:00:00Z"),
                ("2019-11-01T00:00:00Z", "2020-01-01T00:00:00Z")
    */

    var start_date = ee.Date(start);
    var end_date = ee.Date(end);

    // Number of intervals in the given "time_range" based on the specified "agg_interval" period
    var interval_no = ee.Date(end).difference(ee.Date(start), 'day').divide(agg_interval).round();
    var month_check = ee.Algorithms.If(ee.Number(30.4375 / agg_interval).round().gt(0),
                                       ee.Number(30.4375 / agg_interval).round(),
                                       ee.Number(1)); // The number of aggregation intervals within a month

    // Compute the relative date delta (in months) to add to each preceding period to compute the new one
    var rel_delta = ee.Number(end_date.difference(start_date, 'day'))
                    .divide(ee.Number(30.4375).multiply(interval_no)).ceil(); // 30.4375 days = average month length

    // Compute the first time interval end date by adding the relative date delta (in months) to the start date
    end_date = start_date.advance(start_date.advance(rel_delta, 'month')
                                  .difference(start_date, 'day')
                                  .divide(month_check), 'day')
                                  .advance(-1, 'second');

    var time_intervals = ee.List([ee.List([start_date, end_date])]);
    time_intervals = ee.List(ee.List.sequence(1, interval_no.subtract(1)).iterate(function(x,previous){
        start_date = ee.Date(ee.List(ee.List(previous).reverse().get(0)).get(1))
                     .advance(1, 'second'); //end_date of last element
        end_date = start_date
                   .advance(start_date.advance(rel_delta, 'month')
                   .difference(start_date, 'day')
                   .divide(month_check), 'day')
                   .advance(-1, 'second');

        return ee.List(previous).add(ee.List([start_date, end_date]));
    }, time_intervals));

    return time_intervals;
}

exports.makeLegend = function(title, palette, class_names, class_length){
  // Create a legend for the different crop types
  // set position of panel
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '12px 15px'
    }
  });

  // Create legend title
  var legendTitle = ui.Label({
    value: title,
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
      }
  });

  legend.add(legendTitle);

  // Creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
        // Create the label that is actually the colored box.
        var colorBox = ui.Label({
          style: {
            backgroundColor: color,
            // Use padding to give the box height and width.
            padding: '8px',
            fontSize: '12px',
            margin: '0 0 4px 0'
          }
        });

        // Create the label filled with the description text.
        var description = ui.Label({
          value: name,
          style: {margin: '0 0 4px 6px'}
        });

        // return the panel
        return ui.Panel({
          widgets: [colorBox, description],
          layout: ui.Panel.Layout.Flow('horizontal')
        });
  };

  // Add color and and names
  for (var i = 0; i <= class_length; i++) {
    legend.add(makeRow(palette[i], class_names[i]));
    }

  return legend
}

// Function to populate the color palette legends for the app layers
exports.populateLegend = function(legend_name, viz_params, add_char_min, add_char_max, options){

    // Create a legend for the different crop types
    // set position of panel
    var legend = ui.Panel({
      style: {
        position: 'bottom-left',
        padding: '12px 15px'
      }
    });

    // Create legend title
    var legend_title = ui.Label({
      value: legend_name,
      style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 0 0',
      padding: '0',
      //width: '115px'g
      }
      });

    legend.add(legend_title);

    // create the legend image
    var lon = ee.Image.pixelLonLat().select('latitude');
    var gradient = lon.multiply(ee.Number(viz_params.max).subtract(viz_params.min).divide(100)).add(viz_params.min);
    //var viz_palette =  viz_params['palette'].slice(0).reverse();
    //viz_params['palette'] = viz_palette;
    var legend_image = options.legend_image || gradient.visualize(viz_params);
    
    var max_label = ui.Label();
    // create text on top of legend
    var legend_panel_max = ui.Panel({
      widgets: [
      max_label
      //ui.Label(viz_params['max'] + add_char_max)
      ],
      });

    legend.add(legend_panel_max);

    // create thumbnail from the image
    var thumbnail = ui.Thumbnail({
      image: legend_image,
      params: {bbox: '0,0,10,100', dimensions:'10x25'},
      style: {padding: '1px', position: 'bottom-center', fontSize: '18px'}
      });

    legend.add(thumbnail);

    var min_label = ui.Label();
    // create text on top of legend
    var legend_panel_min = ui.Panel({
      widgets: [
      min_label
      ],
      });
      
    ee.Dictionary(viz_params).evaluate(function(params){
      max_label.setValue(params.max + add_char_max)
      min_label.setValue(params.min + add_char_min)
    });
    
    legend.add(legend_panel_min);

    return legend
};
