# Creating your own Wildfire Propagation Database

Code base for the creation of a fire propagation dataset by combining burned area polygons and active fire detection points.
For every burned area, fire propagation is calculated based to active fire detections from remote sensing sensors. The revisiting time of the sensor limits the fire propagation inverval. 
The resulting dataset consists of equally sized raster files, in which each pixel is either classified as burned or not burned. 
Every burned area has at least two distinct raster files, each indicating the burned area at a specific time. With increasing time, more of the final burned area polygon is filled with pixels classified as burned, until every pixel of the area is assined to a fire date.

The following sequence of methods is used to assign a fire date to every pixel of the burned area:
1. The distinct burned area is selected from the burned area dataset. Active fire detections of the same country and the same year are loaded.
2. The dataset of active fire detections is reduced to points that fit a self-defined spatial and temporal range of the burned area polygon.
3. A convex hull is calculated for all active fire detection points of the same aquisition time.
4. Beginning with the youngest date, all convex hulls are iteratively imprinted onto the burned area polygon. This way, pixels with multiple burn dates are overwritten every time a convex hull with another burn date is imprinted, ensuring correct sequencing.
5. Often, not all pixelas are covered by the convex hulls. The remaining pixels are assigned using the k-nearest neighbor algorithm with k=1. This means that these pixels are assigned to the fire date of their geometrically nearest neighbor.
6. Several artifacts can be induced by using the nearest neighbor algorithm. These artifacts are removed by assigning the affected pixels to the next oldest surrounding fire date.
7. In a last step, the calculated fire dates are checked for realistic propagation. If the sequence of propagation is not plausible, the corresponding burned area is discarded. Otherwise, every steps of the fire propagation is saved as an individual raster file.

Most applications of the fire propagation database will involve its combination with different meteorological and/or surface related information. One option to gather such information is Google Earth Engine (GEE), which holds a wide variety of different collections from the earth observation spectrum. Instead of including somewhat arbitrarily chosen datasets directly to the fire propagation database, we decided to provide the necessary code for one meteorological and one remote sensing dataset, i.e. ERA5 and Sentine-2, respectively.
The code enables the download of ERA5 and Sentinel-2 data from GEE for individual burned area polygons, and it can be easily adapted for other GEE-related datasets as well. Providing code instead of real data decreases the physical space of the dataset, while still .  Different examples guide through the creation of the fire propagation database, as well as the download of the corresponding ERA5- and Sentinel-2 data. 


