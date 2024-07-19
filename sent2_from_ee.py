# import necessary libraries
import pandas as pd
import geopandas as geo
import xarray
import rioxarray as rio
from shapely.geometry import box
from osgeo import osr, ogr, gdal
import numpy as np
import rasterio
import datetime
import concurrent
import ee

'''
Module to download Sentinel-2 data from google earth engine.
The pre-processed images are cut to the raster data of the burned areas 
and fit to the burn date.

To ensure cloud free observations, sentinel-2 cloud masking 
from https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless 
was implemented, and cloud-free composite image from successive Sentinel-2 observations are created.

We used the COPERNICUS/S2_HARMONIZED dataset (L1-C) to ensure maximum time coverage, but other datasets can be used as well.
'''

class retrieve_sent2_from_ee:
    
    def __init__(self, burned_area_poly, proj_ee, dataset = 'COPERNICUS/S2_HARMONIZED'):
        self.poly_to_sent = burned_area_poly # geodataframe, in which each row contains a burned area and the associated geometry
        self.dataset = dataset # google earth engine dataset
        self.proj_ee = proj_ee # project name associated with the earth engine account used
        
        
    def authenticate_ee(self):
        '''
        uses your authentification token to activate the python api and connect to your google earth engine account
        
        For further information, read 'https://developers.google.com/earth-engine/guides/auth'
        '''
        proj_ee = self.proj_ee
        
        try:
            ee.Initialize(project=proj_ee, opt_url='https://earthengine-highvolume.googleapis.com')
        except Exception:
            ee.Authenticate()
            ee.Initialize(project=proj_ee, opt_url='https://earthengine-highvolume.googleapis.com')
         
    def get_s2_sr_cld_col(self, aoi, cloud_percentage, start_date, end_date):
        '''
        Functions of the cloud removal algorithm.
        Adapted after 'https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless'
        
        Input:
        aoi = filters collection for the area of interest (ee.geometry). Here, this correponds to the bounding box of a burned area polygon
        cloud_percentage = filters the collection for a define percentage of clouds. 
                        Maximum image cloud cover percent (integer) allowed in image collection.
        start_date = Image collection start date (inclusive; string)
        end_date = Image collection end date (exclusive; string)
        
        Output:
        '''       
        # Import and filter S2 SR.
        s2_sr_col = (ee.ImageCollection(self.dataset) # COPERNICUS/S2_HARMONIZED
            .filterBounds(aoi)
            .filterDate(start_date, end_date)
            .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_percentage)))

        #######################################
        # the original algorithm uses the 'COPERNICUS/S2_SR' dataset, which includes Scene Classification Map.
        # In the original algorithm, this is used to distinguish water bodies from cloud shadows, since both
        # show low reflection values. 
        # COPERNICUS/S2_HARMONIZED does not include such a classification band. In order to identify waterbodies, 
        # # googles dynamic world dataset 'GOOGLE/DYNAMICWORLD/V1' is used and added to the collection.  
        start_date_dynWorld = str(pd.to_datetime(start_date) - datetime.timedelta(days=61))[0:10]
        
        s2_dyn_world = (ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
            .filterBounds(aoi)
            .filterDate(start_date_dynWorld, end_date)
            .mode()
            .select('label'))
        
        def add_new_band(img):
            return img.addBands(s2_dyn_world.rename('LCC'))

        s2_sr_col = s2_sr_col.map(add_new_band)
        #######################################

        # Import and filter s2cloudless.
        s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
            .filterBounds(aoi)
            .filterDate(start_date, end_date))

        # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
        return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
            'primary': s2_sr_col,
            'secondary': s2_cloudless_col,
            'condition': ee.Filter.equals(**{
                'leftField': 'system:index',
                'rightField': 'system:index'
            })
        }))
                
    def add_cloud_bands(self, img):
        # Get s2cloudless image, subset the probability band.
        cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

        cloud_img_thresh = 40
        
        # Condition s2cloudless by the probability threshold value.
        # cloud_img_thresh = Cloud probability (%); values greater than are considered cloud
        is_cloud = cld_prb.gt(cloud_img_thresh).rename('clouds')

        # Add the cloud probability layer and cloud mask as image bands.
        return img.addBands(ee.Image([cld_prb, is_cloud]))

    def add_shadow_bands(self, img):
        # Identify water pixels from the relevant land cover classification band.
        # with 'COPERNICUS/S2_SR', the 'SCL' band can be used.
        # 'COPERNICUS/S2_HARMONIZED' does not have such a band, 
        # and we added the Dynamic World 'LCC' band as a substitution
        
        NIR_thresh = 0.15
        dist_from_cloud = 2
        
        #not_water = img.select('SCL').neq(6)
        not_water = img.select('LCC').neq(0)
        
        ######################

        # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
        # SR_BAND_SCALE is defined according to the employed sentinel-2 dataset, 
        # since Level 1 and Level 2 products can show differences in the band scale
        # the predefined value worked for 'COPERNICUS/S2_HARMONIZED' as well, and, accordingly, we did not change it
        
        # NIR_thresh = Near-infrared reflectance; values less than are considered potential cloud shadow,
        # if they do not belong to a waterbody
        SR_BAND_SCALE = 1e4
        dark_pixels = img.select('B8').lt(NIR_thresh*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

        # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
        shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))

        # Project shadows from clouds for the distance specified by the dist_from_cloud input.
        cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, dist_from_cloud*10)
            .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
            .select('distance')
            .mask()
            .rename('cloud_transform'))

        # Identify the intersection of dark pixels with cloud shadow projection.
        shadows = cld_proj.multiply(dark_pixels).rename('shadows')

        # Add dark pixels, cloud projection, and identified shadows as image bands.
        return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

    def add_cld_shdw_mask(self, img):
        # set buffer around cloud
        cloud_buffer = 100
        # Add cloud component bands.
        img_cloud = self.add_cloud_bands(img)

        # Add cloud shadow component bands.
        img_cloud_shadow = self.add_shadow_bands(img_cloud)

        # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
        is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

        # Remove small cloud-shadow patches and dilate remaining pixels by cloud_buffer input.
        # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
        is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(cloud_buffer*2/20)
            .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
            .rename('cloudmask'))

        # Add the final cloud-shadow mask to the image.
        return img_cloud_shadow.addBands(is_cld_shdw)        

    def apply_cld_shdw_mask(self, img):
        # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
        not_cld_shdw = img.select('cloudmask').Not()

        # Subset reflectance bands and update their masks, return the result.
        return img.select('B.*').updateMask(not_cld_shdw)


    def get_single_images(self, ee_image, index):
        '''
        use earth engines computePixels to get the processed data from server to client.
        This method is significantly faster than former options.
        
        Due to the high resolution of Sentinel-2, the resulting images can be very large.
        To increase speed (and, for large images, even allow the download in the first place), 
        the area of interest is split into multiple smaller areas.
        The smaller areas are then downloaded in parallel, using this function with multiple workers. 
        
        To recreate the original image from the downloaded sub-images, an index is needed to 
        assign each sub-image to its correct position in the original image. This index is therefore part 
        of the input and the output
        
        input: 
        sel_polygon = cloud free composite image to retrieve from the google earth engine server
        index = index to assign the sub-image to its correct position in the large image
        
        output: 
        sel_polygon_npy = numpy array with the defined sentinel-2 bands and the calculated indices 
        index = index to assign the sub-image to its correct position in the large image
        '''    
        try:
            sel_polygon_npy = ee.data.computePixels({
                'expression': ee_image,
                'fileFormat': 'NUMPY_NDARRAY'
            })
        except: 
            sel_polygon_npy = []
        
        return sel_polygon_npy, index
    

    def create_xr_dataset(self, sel_polygon_npy, ee_image, list_of_bands, get_DEM = True, get_LCC = True):
        '''
        function to recreate an xarray from the downloaded numpy arrays.
        The xarray includes the selected bands of sentinel-2.
        If wanted, the DEM and the Dynamic World LCC can be added.
        
        input: 
        sel_polygon_npy = numpy array of the downloaded area
        ee_poly = real google earth engine coordinates of the extracted area. Due to resolution 
        list_of_bands = list with strings of relevant Sentinel-2 bands 
        get_DEM = download the DEM as well. Default is True
        get_LCC = download the Dynamic World LCC as well. Default is True
        
        ouput: 
        data_xr = xarray dataset of the selected area, including all selected Sentinel-2 bands, 
                  as well as DEM and LLC (if True is selected)
        '''
        
        # due to differences in scale and crs, the exact pixels provided by and downloaded from earth engine might differ from the geopandas input
        # to ensure a correct recreation of the image, the outlines of the earth engine image are extracted
        feat_dict = ee_image.getInfo()
        ee_poly = feat_dict['properties']['system:footprint']["coordinates"][0]        
        
        x_new = [ee_poly[0][0], ee_poly[1][0], ee_poly[2][0], ee_poly[3][0], ee_poly[4][0]]
        y_new = [ee_poly[0][1], ee_poly[1][1], ee_poly[2][1], ee_poly[3][1], ee_poly[4][1]]
        
        number_pixel_x = len(sel_polygon_npy[list_of_bands[0]][0])
        number_pixel_y = len(sel_polygon_npy[list_of_bands[0]])
        
        west = min(x_new) + ((max(x_new) - min(x_new))/number_pixel_x)/2
        east = max(x_new) - ((max(x_new) - min(x_new))/number_pixel_x)/2
        south = min(y_new) + ((max(y_new) - min(y_new))/number_pixel_y)/2
        north = max(y_new) - ((max(y_new) - min(y_new))/number_pixel_y)/2

        # use boundaries to create the coordinates of the xarray 
        ncol_coords = np.linspace(start=south, stop=north, num=len(sel_polygon_npy))
        nrow_coords = np.linspace(start=west, stop=east, num=len(sel_polygon_npy[0]))

        # iterate through the list of bands, extract the selected bands and merge them into an xarray dataset
        bands_to_merge = []
        for sel_band in list_of_bands:
            band = sel_polygon_npy[sel_band] 
            band_to_xr = xarray.DataArray(np.flip(band, axis=0), 
                        coords={'y': ncol_coords,'x': nrow_coords}, 
                        dims=["y", "x"])
            band_to_xr.name = sel_band
            band_to_xr = band_to_xr.fillna(np.nan)
            bands_to_merge.append(band_to_xr)
        data_xr = xarray.merge(bands_to_merge)

        # if get_DEM is True, the DEM is extracted and added to the xarray dataset 
        if get_DEM == True:
            dem = sel_polygon_npy["dem30"]
            dem_xr = xarray.DataArray(np.flip(dem, axis=0), 
                        coords={'y': ncol_coords,'x': nrow_coords}, 
                        dims=["y", "x"])
            dem_xr.name = 'dem30'
            dem_xr = dem_xr.fillna(np.nan)

            data_xr =  xarray.merge([data_xr, dem_xr])

        # if get_LCC is True, the Land Cover Classification is extracted and added to the xarray dataset 
        if get_LCC == True:
            lcc = sel_polygon_npy["LCC"]
            lcc_xr = xarray.DataArray(np.flip(lcc, axis=0), 
                        coords={'y': ncol_coords,'x': nrow_coords}, 
                        dims=["y", "x"])
            lcc_xr.name = 'lcc'
            lcc_xr = lcc_xr.fillna(np.nan)
            
            data_xr = xarray.merge([data_xr, lcc_xr])
        
        return data_xr


    def download_sent2(self, sel_polygon, start_date_fire, day_dif_for_composite_img, num_sub_images, out_path, list_of_bands = ["B2", "B3", "B4", "B8", "B11", "B12"], crs = 4326):
        '''
        sel_polygon = geopandas object that contains information about geometry and starting date of a selected burned area (e.g. one row of a geopandas dataframe with equally sized burned area rectangles)
        start_date_fire = starting date of the fire according to the burned area dataset (pandas datetime)
        day_dif_for_composite_img = time difference for which Sentinel-2 images are extracted and the cloud free composit image is prodcued
        num_sub_images = Defines, in how many sub-images the original image is divided. 
                        Depending on the size of the selected burned area polygon,
                        the resulting sentinel-2 images can be too large to transfer from the google earth engine server to the client side.
                        the original image can therefore be split into several sub-images, which are stitched together again after the download.
                        3 images per row and column = 9 sub-images, 4 images per row and column = 16 sub-images
        crs = defines the crs. Default is crs 4326
        out_path = path to save the sentinel-2 .tif-images
        '''
        # all sentinel-2 images between start- and end-date are used to create the cloud-free composit image
        start_date = str(start_date_fire - datetime.timedelta(days = day_dif_for_composite_img))[0:10]
        end_date = str(start_date_fire)[0:10]
        # extract geometry of the burned area
        burned_area_bounds = sel_polygon["geometry"].bounds
        # extract boundary coordinates from the geometry    
        west = burned_area_bounds[0]
        east = burned_area_bounds[2]
        south = burned_area_bounds[1]
        north = burned_area_bounds[3]
        
        new_endpoints_lon = np.linspace(start=west, stop=east, num=(num_sub_images+1))
        new_endpoints_lat = np.linspace(start=south, stop=north, num=(num_sub_images+1))

        new_poly = []
        for long in range(num_sub_images):
            for lati in range(num_sub_images):
                new_poly.append(box(new_endpoints_lon[long], new_endpoints_lat[lati], new_endpoints_lon[long+1], new_endpoints_lat[lati+1]))

        # bounding boxes of the sub-iamges are used to create new geopandas dataframe 
        gdf_poly = geo.GeoDataFrame(new_poly, geometry = new_poly, crs=crs)
        gdf_poly = gdf_poly.drop([0], axis=1)
        
        # iteratively process every sub-image on google earth engines server side
        # processing includes the creation of a cloud-free composite, 
        # as well as adding the DEM and the Dynamic World dataset for the specific region and the specific burn date
        # sub-images are stored in this list
        sub_poly_list = []
        for i in range(len(gdf_poly)):
            selected_poly = gdf_poly.geometry.iloc[i]
            get_coords = list(selected_poly.exterior.coords)
            # define parameters for cloud removal
            aoi = ee.Geometry.Polygon(get_coords)              
            cloud_percentage = 60

            # apply cloud remove 
            s2_sr_cld_col = self.get_s2_sr_cld_col(aoi, cloud_percentage, start_date, end_date)
            s2_sr_median = (s2_sr_cld_col.map(self.add_cld_shdw_mask)
                                    .map(self.apply_cld_shdw_mask)
                                    .median())
            
            # add DEM for specific region
            dem_band = (ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM')
                    .mode()
                    .clip(aoi))
            s2_sr_median = s2_sr_median.addBands(dem_band.rename('dem30'))
            
            # add Land Cover for specific burn date
            # day difference set to 90, to avoid problems with Dynamic Worlds cloud cover
            start_date_dynWorld = str(start_date_fire - datetime.timedelta(days = 90))[0:10]
            s2_dyn_world = (ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
                .filterBounds(aoi)
                .filterDate(start_date_dynWorld, end_date)
                .mode()
                .select('label'))
        
            s2_sr_median = s2_sr_median.addBands(s2_dyn_world.rename('LCC'))
            
            # create the final sub-image after processing
            sub_poly = s2_sr_median.clipToBoundsAndScale(
                geometry=aoi, scale=30
            )
            sub_poly_list.append([sub_poly, i])
        
        # use multiple workers to parallelize the download from server to client side
        # we observed a huge speed up    
        with concurrent.futures.ThreadPoolExecutor(max_workers=100) as executor:
            sub_images = []
            for poly in sub_poly_list:
                sub_images.append(executor.submit(self.get_single_images, ee_image=poly[0], index=poly[1]))
        # save the downloads in a list
        # save the index of each sub-iamge in a different list. 
        # This index is necessary to find the position of the sub-image in the original image        
        list_of_results_init = []
        index_list = []
        for future in concurrent.futures.as_completed(sub_images):
            np_arr, index = future.result()
            index_list.append(index)
            list_of_results_init.append(np_arr)       
        
        # check, if the download was successful
        # if so, recreate original image from sub-images
        if  all(len(x) !=0 for x in list_of_results_init) and all(len(x.dtype) > 2 for x in list_of_results_init):
            
            to_ori_img = [[] for _ in range(num_sub_images)]
            to_shape = 0
            for sub_img in range(0,len(index_list)):
                idx = np.where(np.array(index_list) == sub_img)[0][0]
                if (sub_img == 0) or (sub_img % num_sub_images != 0):
                    to_ori_img[to_shape].append(self.create_xr_dataset(list_of_results_init[idx],sub_poly_list[sub_img][0], list_of_bands = list_of_bands))
                else:
                    to_shape += 1
                    to_ori_img[to_shape].append(self.create_xr_dataset(list_of_results_init[idx],sub_poly_list[sub_img][0], list_of_bands = list_of_bands))
            
            # recreate the original image
            cols_of_ori_img = []
            for sub_col in range(len(to_ori_img[0])):
                cols_of_ori_img.append(xarray.concat(to_ori_img[sub_col], dim='y'))

            ori_img = xarray.concat(cols_of_ori_img, dim='x') 

            # set the crs, convert to rioxarray and save to disk
            ori_img = ori_img.rio.write_crs(f"epsg:{crs}", inplace=True)
            ori_img.rio.to_raster(f"{out_path}/id_{sel_polygon['id']}.tif")
            print(f"{sel_polygon['id']} - successfully saved")
        else:
            # if it did not work, print the problematic id
            print(f"{sel_polygon['id']} - did not work")  
