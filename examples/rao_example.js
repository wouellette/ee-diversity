// ****************************************************************************************************************** //
// ******************************** Rao's Q example for Kalimantan, Indonesia *************************************** //
// ****************************************************************************************************************** //

// Global parameters to define
var START_YEAR = '2018'; // Start year for the Dynamic World Dataset
var END_YEAR = '2022'; // End year for the derivation of temporal indicators and current year indicators
var S2_BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']; // S2 bands to use in the workflow
var S2_COLLECTION = 'COPERNICUS/S2_SR_HARMONIZED'; // S2 collection to use in the workflow
var LS_BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']; // Landsat bands to use in the workflow
var NICFI_BANDS = ['B', 'G', 'R', 'N'];
var SURFACE_WATER = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select('max_extent').eq(1); // Surface Water Mask
var CLD_PRB_THRESH = 20; // Cloud probability threshold to mask clouds. 40% is the default value of s2cloudless
var CLD_FILTER = 40; // Threshold on sentinel-2 Metadata field determining whether cloud pixel percentage in image
var NIR_DRK_THRESH = 0.15; // A threshold that determines when to consider a dark area a cloud shadow or not
var CLD_PRJ_DIST = 10; // The distance (in no of pixels) in which to search from detected cloud to find cloud shadows
var CLD_BUFFER = 200; // The cloud buffer (in meters) to use around detected cloud pixels to mask additionally
var SNW_THRESH = 0.4; // The snow buffer (in meters) to use around detected cloud pixels to mask additionally
var MASK_RES = 60; // resolution at which to generate and apply the cloud/shadow mask. 60m instead of 10m to speed up
var AOI = //ee.FeatureCollection("projects/sat-io/open-datasets/geoboundaries/CGAZ_ADM2")
          //.filter(ee.Filter.eq('shapeName', 'Ogan Komering Ulu Selatan'));
          ee.FeatureCollection(ee.Feature(ee.Geometry.Rectangle([104.0810472285037,-4.62650507909588,
                                                                 104.10967176311064,-4.609052455388382]), null));

// Rao's Q specific parameters
var WINDOW_SIZE = 9; // Default window size when applied to Sentinel-2 resolution data
var ALPHA = 2; // Quadratic alpha parameter

// Dependencies Import
var palettes = require('users/gena/packages:palettes');
var utils = require('users/soilwatchtech/biodiversityApp:utils.js');
var diversity = require('users/soilwatchtech/biodiversityApp:rao.js');

var dynamicWorld_class_names = [
  'built', 'crops', 'trees', 'shrub_and_scrub', 'grass',
  'water', 'flooded_vegetation', 'bare', 'snow_and_ice'];
  // The corresponding color hex keys for the land cover classes
var dynamicWorld_palette =
  [
  '#c4281b', // 1.urban
  '#e49634', // 2.croplands
  '#397e48', // 3. forest
  '#DFC35A', // 4. shrub & scrub
  '#88af52', // 5. grass
  '#429ae4', // 6. water
  '#7c85c9', // 7. wetlands
  '#a6a6a6', // 8. bare surface
  '#B39FE1' // 9. ice and snow
  ];

var dynamicWorld_legend = utils.makeLegend('Dynamic World', dynamicWorld_palette, dynamicWorld_class_names, 8);
Map.add(dynamicWorld_legend);

// Load latest available land cover datasets (World Cover, Dynamic World)
var landcover_data = utils.landCoverDatasets(END_YEAR, AOI);
var land_cover = landcover_data[0];
var tree_cover = landcover_data[1];
var wetland_cover = landcover_data[2];
var dynamic_world = landcover_data[3];
var proba_hillshade = landcover_data[4];
var glad_change = landcover_data[5];
var forest_change = landcover_data[6];

// Plot Dynamic World data aggregated for 2022
Map.centerObject(AOI);
Map.setOptions("SATELLITE");
Map.addLayer(AOI, {}, 'AOI');
Map.addLayer(dynamic_world.updateMask(proba_hillshade.divide(255).gt(0.65)).clip(AOI.geometry())
             .visualize({min: 1, max: dynamicWorld_palette.length, palette: dynamicWorld_palette})
             .divide(255).multiply(proba_hillshade.divide(255)),
             {min: 0, max: 0.65},
             'Dynamic World Land Cover 2022');

// ******************** spatio-parametric Rao's Q for single date applied to Sentinel-2 NDVI ************************ //

var s2_col = utils.S2CloudMasked(S2_COLLECTION, S2_BANDS, START_YEAR, END_YEAR, SURFACE_WATER,
                                        CLD_FILTER, NIR_DRK_THRESH, CLD_PRJ_DIST, CLD_PRB_THRESH, CLD_BUFFER,
                                        SNW_THRESH, MASK_RES, AOI)
                                        .map(utils.applyBRDF('s2'));

// Compute the max ndvi (95th percentile) for the current year
var ndvi_median_s2 = ee.ImageCollection(ee.List.sequence(ee.Number.parse(START_YEAR), ee.Number.parse(END_YEAR))
                     .map(function(year){
                       year = ee.String(ee.Number(year).toInt16());
                       var ndvi_median = s2_col.filterDate(year.cat('-01-01'), year.cat('-12-31'))
                        .map(function(img){return img.divide(10000).addBands(img.normalizedDifference(['B8', 'B4'])
                                                     .rename('NDVI')
                                                     .copyProperties(img, ['system:time_start']))
                        })
                        .reduce(ee.Reducer.percentile([50]), 4)
                        .rename(S2_BANDS.concat(['NDVI'])).toFloat()
                        .set('system:time_start', ee.Date(year.cat('-07-01')).millis());

                        return ndvi_median
                     }));

// Get the latest year of Sentinel-2 NDVI data
var ndvi_latest_s2 = ndvi_median_s2.sort('system:time_start', false).first();

// Get the most recent Rao's Q entropy calculation
var rao_s2_ndvi_col = ndvi_median_s2.select('NDVI') // Embed in image coll as RaoQ is wrapper function
                    .map(diversity.raoQ(AOI.geometry(), {window_size: WINDOW_SIZE, alpha: ALPHA}))
                    .map(utils.createTimeBand);

var rao_s2_ndvi = rao_s2_ndvi_col.sort('system:time_start', false).first();

// Get the most recent Rao's Q entropy calculation
var rao_s2_multi_col = ndvi_median_s2.select(S2_BANDS) // Embed in image coll as RaoQ is wrapper function
                    .map(diversity.raoQ(AOI.geometry(), {window_size: WINDOW_SIZE, alpha: ALPHA}))
                    .map(utils.createTimeBand);

var rao_s2_multi = rao_s2_multi_col
                   .sort('system:time_start', false).first()
                   .reduce(ee.Reducer.sum()); // Additive property of Rao's Q enables summation of the bands

// Load pre-computed results for faster loading in code editor
rao_s2_ndvi = ee.Image('users/soilwatchtech/Indonesia/rao_s2_ndvi_2022');
//rao_s2_multi = ee.Image('users/soilwatchtech/Indonesia/rao_s2_multi_2022');

// Load color palettes from GEE community palettes
var rao_palette = palettes.matplotlib.viridis[7].slice(0);
var rao_change_palette = palettes.colorbrewer.RdYlGn[11].slice(0);
var rao_change_palette_rev = palettes.colorbrewer.RdYlGn[9].slice(0).reverse();

var raoLatest_legend = utils.populateLegend('Beta diversity 2022',
                                  {min: 0, max: 1, palette: rao_palette},
                                  " (low)", " (high)", {});
Map.add(raoLatest_legend);

var ndvi_legend = utils.populateLegend('NDVI/CV/trend',
                                  {min: 0, max: 1, palette: rao_change_palette},
                                  " -", " +", {});
Map.add(ndvi_legend);

Map.addLayer(ndvi_latest_s2.select(['B4', 'B3', 'B2']).clip(AOI.geometry()),
             {min: 0, max: 0.1},
             "Median Sentinel-2 RGB 2022");

Map.addLayer(ndvi_latest_s2.select('NDVI').clip(AOI.geometry()),
             {min: 0.4, max: 0.9, palette: rao_change_palette},
             "Median Sentinel-2 NDVI 2022");

Map.addLayer(rao_s2_ndvi.updateMask(rao_s2_ndvi.lte(0.7) // Mask above 0.7, corresponds to forest/non-forest interface
                                    .and(ndvi_latest_s2.select('NDVI').gt(0.8)) // High NDVI values (i.e. closed canopy)
                                    .and(dynamic_world.eq(3)) // Dynamic world forest class
                                    .and(proba_hillshade.divide(255).gt(0.65))) // Tree proba > 65% as per Dynamic World
             .reproject({scale:10, crs:'EPSG:32748'})
             .clip(AOI.geometry()),
             {min: 0, max: 1, palette: rao_palette},
             "Rao's Q (alpha=2, window_size=9) Sentinel-2 NDVI");


Map.addLayer(rao_s2_multi.updateMask(rao_s2_multi.lte(3.5)
                                       .and(ndvi_latest_s2.select('NDVI').gt(0.8))
                                       .and(dynamic_world.eq(3))
                                       .and(proba_hillshade.divide(255).gt(0.65)))
             .reproject({scale:10, crs:'EPSG:32748'})
             .clip(AOI.geometry()),
             {min: 0, max: 4, palette: rao_palette},
             "Rao's Q (alpha=2, window_size=9) Sentinel-2 Multi-band");

// ***************** spatio-parametric Rao's Q for single date applied to a canopy height model ********************* //

// Load Lang et al., 2021's canopy height model for the year 2020
var canopy_height = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1').rename('canopy_height');

// Calculate Rao's Q on Lang et al., 2021's canopy height model
var rao_ch = ee.ImageCollection([canopy_height.clamp(0, 43) // Clamping and scaling canopy height to [0,1] range
                                              .unitScale(0, 43)]) // Embed in imageColl as RaoQ is wrapper
             .map(diversity.raoQ(AOI.geometry(), {window_size: WINDOW_SIZE, alpha: ALPHA}))
             .first();

// Load pre-computed results for faster loading in code editor
rao_ch = ee.Image('users/soilwatchtech/Indonesia/rao_ch_2020');

Map.addLayer(canopy_height.clip(AOI.geometry()),
             {min: 0, max: 43, palette: rao_palette},
             "Lang et al., 2021 10m canopy height model");

Map.addLayer(rao_ch.updateMask(ndvi_latest_s2.select('NDVI').gt(0.8)
                               .and(dynamic_world.eq(3))
                               .and(proba_hillshade.divide(255).gt(0.65)))
             .reproject({scale:10, crs:'EPSG:32748'})
             .clip(AOI.geometry()),
             {min: 0, max: 2, palette: rao_palette},
             "Rao's Q (alpha=2, window_size=9) Lang et al., 2021 Canopy Height Model");

// Summing rescaled Rao's Q outputs for Sentinel-2 NDVI and Lang et al.,2021 Canopy Height as generalized Additive Model
Map.addLayer(rao_ch.clamp(0, 2).unitScale(0, 2).add(rao_s2_ndvi.clamp(0, 1))
             .updateMask(rao_s2_ndvi.lte(0.7)
                         .and(ndvi_latest_s2.select('NDVI').gt(0.8))
                         .and(dynamic_world.eq(3))
                         .and(proba_hillshade.divide(255).gt(0.65)))
             .reproject({scale:10, crs:'EPSG:32748'})
             .clip(AOI.geometry()),
             {min: 0, max: 2, palette: rao_palette},
             "Generalized Additive Model (Rao's Q S2 NDVI + Rao's Q Canopy Height)")

// ********************* spatio-parametric Rao's Q applied to a time series of NICFI Basemaps *********************** //

// Load NICFI Basemaps for Asia (Indonesia)
var nicfi_col = ee.ImageCollection("projects/planet-nicfi/assets/basemaps/asia")
        .select(NICFI_BANDS)
        .filterBounds(AOI.geometry())
        .map(function(img){return img.divide(10000).addBands(img.normalizedDifference(['N', 'R']).rename('NDVI'))
                                  .copyProperties(img, ['system:time_start']);
        });

// Extract the median from the NICFI basemaps for a given calendar year
var nicfi_col_ndvi = ee.ImageCollection(ee.List.sequence(ee.Number.parse(START_YEAR), ee.Number.parse(END_YEAR))
                     .map(function(year){
                       year = ee.String(ee.Number(year).toInt16());
                        var ndvi_max_nicfi = nicfi_col.select('NDVI').filterDate(year.cat('-01-01'), year.cat('-12-31'))
                                             .reduce(ee.Reducer.percentile([50]), 16)
                                             .rename('NDVI').toFloat()
                                             .set('system:time_start', ee.Date(year.cat('-07-01')).millis());
                        return ndvi_max_nicfi;
                     }));

// Get the time series of Rao's Q entropy calculation from the NICFI NDVI time series
var rao_ndvi_nicfi = nicfi_col_ndvi
                     .map(diversity.raoQ(AOI.geometry(), {window_size: WINDOW_SIZE*2, alpha: ALPHA}))
                     .map(utils.createTimeBand);

// Get the time series of Rao's Q entropy calculation from the NICFI spectral bands time series
// Extract the corresponding spectral bands to the max NDVI from the NICFI basemaps for a given calendar year
/*var nicfi_col_multi = ee.ImageCollection(ee.List.sequence(ee.Number.parse(START_YEAR), ee.Number.parse(END_YEAR))
                     .map(function(year){
                       year = ee.String(ee.Number(year).toInt16());
                        var ndvi_max_nicfi = nicfi_col.filterDate(year.cat('-01-01'), year.cat('-12-31'))
                                             .qualityMosaic('NDVI')
                                             .select(NICFI_BANDS).toFloat()
                                             .set('system:time_start', ee.Date(year.cat('-07-01')).millis());
                        return ndvi_max_nicfi;
                     }));

var rao_multi_nicfi = nicfi_col_multi
                      .map(diversity.raoQ(AOI.geometry(), {window_size: WINDOW_SIZE, alpha: ALPHA}))
                      .map(function(img){
                        return img.reduce(ee.Reducer.sum())
                                  .copyProperties(img, ['system:time_start'])})
                      .map(utils.createTimeBand);*/

// Derive the Rao Coefficient of Variation
var rao_latest = rao_ndvi_nicfi.sort('system:time_start', false).first().select('NDVI');
var rao_mean = rao_s2_ndvi_col.select('NDVI').reduce(ee.Reducer.mean(), 16);
var rao_std = rao_s2_ndvi_col.select('NDVI').reduce(ee.Reducer.stdDev(), 16);
var rao_cv = rao_std.divide(rao_mean).rename('rao_ndvi_cv');

var rao_ols = utils.trendTS(rao_s2_ndvi_col.select('NDVI'), 'system:time_start', 'NDVI');

// Load pre-computed results for faster loading in code editor
rao_latest = ee.Image('users/soilwatchtech/Indonesia/rao_nicfi_ndvi_2022');
rao_ols = ee.Image('users/soilwatchtech/Indonesia/rao_s2_trend_2022');
rao_cv = ee.Image('users/soilwatchtech/Indonesia/rao_s2_cv_2022');

Map.addLayer(nicfi_col.sort('system:time_start', false).first().select(['R', 'G', 'B']).clip(AOI.geometry()),
             {min: 0, max: 0.1},
             'Median NICFI RGB 2022');

Map.addLayer(nicfi_col_ndvi.sort('system:time_start', false).first().clip(AOI.geometry()),
             {palette: rao_change_palette, min: 0.4, max: 0.9},
             'Median NICFI NDVI 2022');

Map.addLayer(rao_latest
             .updateMask(rao_latest.lte(0.7)
                         .and(nicfi_col_ndvi.sort('system:time_start', false).first().gt(0.8))
                         .and(dynamic_world.eq(3))
                         .and(proba_hillshade.divide(255).gt(0.65)))
             .clip(AOI.geometry()),
             {palette: rao_palette, min: 0, max: 1},
             "Rao's Q NICFI NDVI 2022 (alpha=2, window_size=18)");

Map.addLayer(rao_ols.select('scale')
                    .updateMask(rao_s2_ndvi.lte(0.7)
                         .and(ndvi_latest_s2.select('NDVI').gt(0.8))
                         .and(dynamic_world.eq(3))
                         .and(proba_hillshade.divide(255).gt(0.65)))
                    .clip(AOI.geometry()),
             {palette: rao_change_palette, min: -0.1, max: 0.1},
             "Rao's Q S2 NDVI slope of change 2018-2022 (alpha=2, window_size=18)");

Map.addLayer(rao_cv.updateMask(rao_s2_ndvi.lte(0.7)
                         .and(ndvi_latest_s2.select('NDVI').gt(0.8))
                         .and(dynamic_world.eq(3))
                         .and(proba_hillshade.divide(255).gt(0.65)))
                    .clip(AOI.geometry()),
             {palette: rao_change_palette_rev, min: 0, max: 1},
             "Rao's Q S2 NDVI Coefficient of Variation 2018-2022 (alpha=2, window_size=18)");

/*
// Define arguments for animation function parameters.
var video_args = {
  'dimensions': 2400,
  'region': AOI.geometry(),
  'framesPerSecond': 1,
  'crs': 'EPSG:32748',
  'min': 0,
  'max': 1,
  'palette': rao_palette
}

// Function to export GIF
var generateGIF = function(img, band_name, args){
  print(img.select(band_name).getVideoThumbURL(args));
}

// Export GIF of the Rao's Q timelapse
print("Rao's Q NICFI timelapse:");
generateGIF(rao_s2_ndvi, 'NDVI', video_args);
*/

// Export generated Rao's Q results. This is recommended to get the full extent of the data generated and saved,
// so it can be explored and consulted seamlessly.
Export.image.toAsset({
  image: rao_s2_multi.clip(AOI.geometry()),
  description:'rao_s2_multi_2022',
  assetId: 'users/soilwatchtech/Indonesia/rao_s2_multi_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 10,
  maxPixels:1e13
});

Export.image.toAsset({
  image: rao_s2_ndvi.clip(AOI.geometry()),
  description:'rao_s2_ndvi_2022',
  assetId: 'users/soilwatchtech/Indonesia/rao_s2_ndvi_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 10,
  maxPixels:1e13
});

Export.image.toAsset({
  image: rao_ch.clip(AOI.geometry()),
  description:'rao_ch_2020',
  assetId: 'users/soilwatchtech/Indonesia/rao_ch_2020',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 10,
  maxPixels:1e13
});

Export.image.toAsset({
  image: rao_ndvi_nicfi.sort('system:time_start', false).first().select('NDVI').clip(AOI.geometry()),
  description:'rao_nicfi_ndvi_2022',
  assetId: 'users/soilwatchtech/Indonesia/rao_nicfi_ndvi_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

Export.image.toAsset({
  image: rao_cv.clip(AOI.geometry()),
  description:'rao_nicfi_cv_2022',
  assetId: 'users/soilwatchtech/Indonesia/rao_s2_cv_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

Export.image.toAsset({
  image: rao_ols.clip(AOI.geometry()),
  description:'rao_nicfi_trend_2022',
  assetId: 'users/soilwatchtech/Indonesia/rao_s2_trend_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});
