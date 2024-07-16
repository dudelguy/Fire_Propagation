# Creating your own Wildfire Propagation Database

Code base for the creation of a fire propagation dataset by combining burned area polygons and active fire detection points.
For every burned area, fire propagation is calculated based to active fire detections from remote sensing sensors. The revisiting time of the sensor limits the fire propagation inverval. 
The resulting dataset consists of equally sized raster files, in which each pixel is either classified as burned or not burned. 
Every burned area has at least two distinct raster files, each indicating the burned area at a specific time. With increasing time, more of the final burned area polygon is filled with pixels classified as burned, until every pixel of the area is assined to a fire date.

The following sequence of methods is used to assign a fire date to every pixel of the burned area:
1. The distinct burned area is selected from the burned area dataset. Active fire detections of the same country and the same year are loaded.

![step_1](https://github.com/user-attachments/assets/63540813-191f-4f57-aa86-a5a2a735e31e)


2. The dataset of active fire detections is reduced to points that fit a self-defined spatial and temporal range of the burned area polygon.

![step_2](https://github.com/user-attachments/assets/9f58e4da-3b4f-44c7-8aa3-dfb3df9115a7)

3. A convex hull is calculated for all active fire detection points of the same aquisition time.

![step_3](https://github.com/user-attachments/assets/98756863-11c5-411f-a918-49182d4173a0)

4. Beginning with the youngest date, all convex hulls are iteratively imprinted onto the burned area polygon. This way, pixels with multiple burn dates are overwritten every time a convex hull with another burn date is imprinted, ensuring correct sequencing.

![conv_hull_imprint](https://github.com/user-attachments/assets/2a045bd7-144d-4838-a16e-13b6dee3ed05)

5. Often, not all pixelas are covered by the convex hulls. The remaining pixels are assigned using the k-nearest neighbor algorithm with k=1. This means that these pixels are assigned to the fire date of their geometrically nearest neighbor.

![step_5](https://github.com/user-attachments/assets/5245e337-4175-4c66-bfcb-673f73c4b97e)

6. Several artifacts can be induced by using the nearest neighbor algorithm. These artifacts are removed by assigning the affected pixels to the next oldest surrounding fire date.

![step_6](https://github.com/user-attachments/assets/4c152b9d-7786-441f-b7d9-8ce78156ccd1)

7. In a last step, the calculated fire dates are checked for realistic propagation. If the sequence of propagation is not plausible, the corresponding burned area is discarded. Otherwise, every steps of the fire propagation is saved as an individual raster file.

![final_propagation](https://github.com/user-attachments/assets/23264620-68b8-4905-b4ab-1a58927965a4)



Most applications of the fire propagation database will involve its combination with different meteorological and/or surface related information. One option to gather such information is Google Earth Engine (GEE), which holds a wide variety of different collections from the earth observation spectrum. Instead of including somewhat arbitrarily chosen datasets directly to the fire propagation database, we decided to provide the necessary code for one meteorological and one remote sensing dataset, i.e. ERA5 and Sentine-2, respectively.
The code enables the download of ERA5 and Sentinel-2 data from GEE for individual burned area polygons, and it can be easily adapted for other GEE-related datasets as well. Providing code instead of real data decreases the physical space of the dataset, while still .  Different examples guide through the creation of the fire propagation database, as well as the download of the corresponding ERA5- and Sentinel-2 data. 


