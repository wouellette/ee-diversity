// ****************************************************************************************************************** //
// ****************************** Function to calculate spatio-parametric Rao's Q *********************************** //
// ****************************************************************************************************************** //

/**
 * Wrapper function to calculate the spatio-temporal parametric rao's Q as per Rocchini et al., 2021,
   as a proxy for Plant Species Diversity (Beta Diversity)
 * @param {ImageCollection} collection: The input satellite Collection (single or multi-band),
                                        to which a mapping function should be applied
 * @param {Number} window_size: The window size to use for the spatial computation
 * @param {Number} alpha: The alpha parameter of the parametric Rao's Q
 * @returns {ImageCollection}: A time series for which the Rao's Q was calculated at each timestamp and for each band.
                               The resulting bands can be summed a posteriori to generate a single Rao's Q index
                               due to its additive property
 * @ignore
*/
exports.raoQ = function(options){

    var window_size = options.window_size || 9; // Typical window size of 9 is advised by authors for Sentinel-2
    var alpha = options.alpha || 2; // 2 is the quadratic form of Rao's Q

    // Apply wrapper function to allow applying Rao's Q to time series data
    var wrap = function(img){
      // Iterate over the input bands
      var rao_multi = ee.ImageCollection(img.bandNames().map(function(band_name){
          // Extract the moving window for each center pixel
          var subsets = img.select([band_name])
                           .neighborhoodToArray(ee.Kernel.square(window_size/2));
          // Number of pixels in the moving window
          var N = ee.Number(window_size).pow(2);

          // Compute the terms of Rao's Q to sum up (distance matrix and generalized mean)
          // Pairwise pixel distances are computed in both dimensions of the moving window
          var rao_q = ee.Image(ee.List.sequence(1, window_size-2).iterate(function(i, img) {
            return ee.Image(img).add(ee.List.sequence(ee.Number(i).add(1), window_size-1).iterate(function(j, img) {
              return ee.Image(img).add(subsets
                                      .subtract(subsets.arraySlice(0, ee.Number(i).int(), ee.Number(i).add(1).int())
                                                       .arraySlice(1, ee.Number(j).int(), ee.Number(j).add(1).int())
                                                       .arrayProject([0])
                                                       .arrayFlatten([['community']])
                                               ).abs().pow(alpha)
                                      .multiply(ee.Image(1).divide(N.pow(2)))
                                      .pow(ee.Number(1).divide(alpha))
                                      );
            }, subsets.subtract(subsets.arraySlice(0, ee.Number(i).int(), ee.Number(i).add(1).int())
                                       .arraySlice(1, 1, 2)
                                       .arrayProject([0])
                                       .arrayFlatten([['community']])
                               ).abs().pow(alpha)
                               .multiply(ee.Image(1).divide(N.pow(2))
                               .pow(ee.Number(1).divide(alpha))
                               )));
          }, subsets.subtract(subsets.arraySlice(0, 0, 1)
                                     .arraySlice(1, 1, 2)
                                     .arrayProject([0])
                                     .arrayFlatten([['community']])
                              ).abs().pow(alpha)
                              .multiply(ee.Image(1).divide(N.pow(2)))
                              .pow(ee.Number(1).divide(alpha))
          ));

          return rao_q.arrayReduce({reducer:ee.Reducer.sum().unweighted(), axes:[0, 1]}) // Sum all moving window terms
                      .arrayProject([0])
                      .arrayFlatten([[ee.String(band_name).cat('_rao')]]) // Convert array to image with band name
                      .divide(2) // Divide by two because of distance matrix symmetry
                      .copyProperties(img, ['system:time_start']); // Carry over system:time_start from original image
      }));

      return rao_multi.toBands().rename(img.bandNames()).copyProperties(rao_multi.first(), ['system:time_start']);
    }

    return wrap
};
