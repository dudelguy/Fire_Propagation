# Fire_Propagation

Code base for the creation of a fire propagation dataset by combining burned area shapefiles and active fire detections.
For every burned area, fire propagation is calculated based to active fire detections from remote sensing sensors. The revisiting time of the sensor limits the fire propagation inverval. 
The resulting dataset consists of equally sized raster files, in which each pixel is either classified as burned or not burned. 
Every burned area has at least two distinct raster files, each indicating the burned area at a specific time. With increasing time, more of the final burned area polygon is filled with pixels classified as burned, until every pixel of the area is assined to a fire date.

The following sequence of methods is used to assign a fire date to every pixel of the burned area:
1. The distinct burned area is selected.
2. All active fire detection points in a predefined spatial and temporal range of the burned area polygon are extracted.
3. A convex hull is calculated for all points of the same aquisition time.
4. Beginning with the youngest date, all convex hulls are iteratively imprinted onto the burned area polygon. This way, pixels with multiple burn dates are overwritten every time a convex hull with an oder burn date is imprinted, ensuring correct sequencing.
5. The remaining pixels that are assigned using the k nearest neighbor algorithm with k=1. This means that these pixels are assigned to the fire date of their geometrically nearest neighbor.
6. Several artifacts can be induced by using the nearest neighbor algorithm. These artifacts are removed by assigning the affected pixels to the next oldest surrounding fire date.
7. In a last step, the calculated fire dates are checked for realistic propgation. If the sequence of propagation is not plausible, the corresponding burned area is discarded and no raster files are saved. Otherwise, every steps of the fire propagation is saved as an individual raster file.

Using the fire propagation database involves a combination with different meteorological and/or surface related informations. Google Earth Engine (GEE) holds a wide variety of different collections, including ERA5 reanalysis data and Sentinel-2 images of different processing stages. Instead of including this data directly to th fire propagation database, the necessary code to download ERA5 and Sentine-2 data from GEE for individual burned area polygons is provided. This prevents unecessary fire propagation images. Different examples guide through the creation of the fire propagation database, as well as the download of the corresponding ERA5- and Sentinel-2 data. 
