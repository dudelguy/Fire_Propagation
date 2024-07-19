# Creating your own Wildfire Propagation Database

This repository contains the code base for the creation of a fire propagation dataset by combining burned area polygons and active fire detection points.
For every burned area, fire propagation is calculated based to active fire detections from remote sensing sensors. The revisiting time of the sensor limits the fire propagation interval. 
The resulting dataset consists of equally sized raster files, in which each pixel is either classified as burned or not burned. 
Every burned area is associated with at least two distinct raster files, each indicating the burned area at a specific time. With increasing time, more of the final burned area polygon is filled with pixels classified as burned, until every pixel of the area is assigned to a fire date.

The following sequence of methods is used to assign a fire date to every pixel of the burned area:

### 1. Load the data
The distinct burned area is selected from the burned area dataset. Active fire detections of the same country and the same year are loaded.

<p align="center">
<img src="https://github.com/user-attachments/assets/88f6760f-f0c2-47ca-a6de-f8d0c874a300" width=75% height=75%>
</p>

### 2. Reduce active fire detections 
The dataset of active fire detections is reduced to points that fit a self-defined spatial and temporal range of the burned area polygon.

<p align="center">
<img src="https://github.com/user-attachments/assets/4db1bc37-053d-4472-9df9-6055a6050dc4" width=75% height=75%>
</p>

### 3. Calculate convex hull
A convex hull is calculated for all active fire detection points of the same acquisition time. According to the active fire detection algorithm developed by [Schroeder and Gigilio 2016](https://viirsland.gsfc.nasa.gov/PDF/VIIRS_activefire_375m_ATBD.pdf), thermal anomalies are based on temperature differences of the investigated pixel to a defined number of surrounding pixels. If larger areas are on fire, all surrounding pixels show large temperature values as well, and, as a result, no anomaly will be detected and the pixel will not be classified as an active fire pixel. In such cases, only the edges of a fire area will be detected as active, although the area in between potentially burns, too. To account for these situations, the convex hull was calculated and all affected pixel inside of this hull were assumed to be burning.  

<p align="center">
<img src="https://github.com/user-attachments/assets/57a57c88-137b-4ff2-b7e3-3c684f93d332" width=65% height=65%>
</p>

### 4. Fill pixels according to convex hulls
Beginning with the youngest date, all convex hulls are iteratively imprinted onto the burned area polygon. This way, pixels with multiple burn dates are overwritten every time a convex hull with an older burn date is imprinted, ensuring correct sequencing. 

<p align="center">
<img src="https://github.com/user-attachments/assets/b11be3a4-c45b-4da3-a871-8a9a8268a645" width=50% height=50%>
</p>

### 5. Fill remaining pixels with k-nearest neighbor
Often, not all pixelas are covered by the convex hulls. The remaining pixels are assigned using the k-nearest neighbor algorithm with k=1. This means that these pixels are assigned to the fire date of their geometrically nearest neighbor.

<p align="center">
<img src="https://github.com/user-attachments/assets/12b49eeb-a09e-46bd-8c7b-ba868ef98ed1" width=50% height=50%>
</p>

### 6. Remove artifacts
Several artifacts can be appear due to the use of the nearest neighbor algorithm. These artifacts are removed by assigning the affected pixels to the next oldest surrounding fire date.

<p align="center">
<img src="https://github.com/user-attachments/assets/a7a88e2b-c2c0-4e04-863b-ac63463b4c90" width=50% height=50%>
</p>

### 7. Check fire propagation
In a last step, the calculated fire dates are checked for realistic propagation. If the sequence of propagation is not plausible, the corresponding burned area is discarded. Otherwise, every steps of the fire propagation is saved as an individual raster file.

<p align="center">
<img src="https://github.com/user-attachments/assets/e6aee86e-eb86-4d2f-8bde-211b18c479e6">
</p>

# Matching the burned areas with surface information and meteorological data

Most applications of the fire propagation database will involve its combination with different meteorological and/or surface related information. One option to gather such information is Google Earth Engine (GEE), which holds a wide variety of different collections from the earth observation spectrum. Instead of including somewhat arbitrarily chosen datasets directly to the fire propagation database, we decided to provide the necessary code for one meteorological and one remote sensing dataset, i.e. ERA5 and Sentine-2, respectively.
The code enables the download of ERA5 and Sentinel-2 data from GEE for individual burned area polygons, and it can be easily adapted for other GEE-related datasets as well. Providing code instead of real data decreases the physical space of the dataset, while still .  Different examples guide through the creation of the fire propagation database, as well as the download of the corresponding ERA5- and Sentinel-2 data. 


