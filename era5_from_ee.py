
# import necessary libraries
import pandas as pd
import geopandas as geo
import xarray
import rioxarray as rio
from rasterio import Affine
from osgeo import osr, ogr, gdal
import os
import numpy as np
import rasterio
import datetime
import concurrent
import ee

'''
Module to download ERA5-reanalysis data from google earth engine.
The reanalysis results are cut to the size of the burned areas and fit to the burn date.

ERA5 offers hourly updates with low spatial resolution (9km). 
Meteorological data is tehrefore available for every hour of the wildfire.
'''

class retrieve_era5_from_ee:
    
    def __init__(self, burned_area_poly, proj_ee, dataset = 'ECMWF/ERA5_LAND/HOURLY'):
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
            
    def wind_speed_era5(self, u_speed, v_speed):
        '''
        ERA5 returns wind data as its u- and v-components.
        This function can be used to calculate the according wind speed
        
        input:
        u_speed = u-component of the wind
        v_speed = v-component of the wind
        
        output: 
        wind speed along the direction of the wind 
        '''
        return np.sqrt((u_speed**2)+(v_speed**2))
    
    
    def wind_dir_era5(self, u_speed, v_speed):
        '''
        ERA5 returns wind data as its u- and v-components.
        This function can be used to calculate the wind direction (in degrees FROM which the wind is coming from)
        
        input:
        u_speed = u-component of the wind
        v_speed = v-component of the wind
        
        output: 
        direction, from where the wind is coming from (degrees from 0-360)
        '''
        return (180+((180/np.pi)*np.arctan2(u_speed,v_speed)))%360
    
    
    def calculate_relative_humidity(self, T, Dp):
        '''
        ERA5 has no parameter for the relative humidity.
        This function can be used to calculate the relative humidity from the ERA5 temperature (in K) and dew temperature values
        
        input:
        T = Temperature in K
        Dp = Dew temperature
        
        output: relative humidity
        '''
        b = 17.625
        c = 243.04

        numerator = np.exp((b * (Dp-273.15)) / (c + (Dp-273.15)))
        denominator = np.exp((b * (T-273.15)) / (c + (T-273.15)))

        RH = 100 * (numerator / denominator)
        return RH
    
    
    def get_single_images(self, ee_image, index):
        '''
        use earth engines computePixels to get the processed data from server to client.
        This method is significantly faster than former options.
        
        input: 
        sel_polygon = image to retrieve from the google earth engine server
        index = index to assign the sub-image to its correct position in the large image
        
        output: 
        sel_polygon_npy = numpy array with the defined ERA5 component
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


    def create_xr_dataset(self, sel_polygon_npy, ee_image, list_of_bands, range_hours):
        '''
        function to recreate an xarray from the downloaded numpy arrays.
        The xarray includes the selected bands of sentinel-2.
        If wanted, the DEM and the Dynamic World LCC can be added.
        
        input: 
        sel_polygon_npy = numpy array of the downloaded area
        ee_poly = real google earth engine coordinates of the extracted area. Due to resolution 
        list_of_bands = list with strings of relevant Sentinel-2 bands 
        range_hours = range of hours for which to download the ERA5 data
        
        ouput: 
        data_xr = xarray dataset of the selected area, including all selected Sentinel-2 bands, 
                  as well as DEM and LLC (if True is selected)
        '''
        
        '''
        # get the boundaries of the downloaded area from the original geopandas dataframe
        west = selected_poly.bounds[0]
        east = selected_poly.bounds[2]
        south = selected_poly.bounds[1]
        north = selected_poly.bounds[3]
        '''
        # due to differences in scale and crs, the exact pixels provided by and downloaded from earth engine might differ from the geopandas input
        # to ensure a correct recreation of the image, the outlines of the earth engine image are extracted
        feat_dict = ee_image.getInfo()
        ee_poly = feat_dict['properties']['system:footprint']["coordinates"][0]        
        
        # x- and y-coordinates of every pixel are created using the earth engine information
        x_new = [ee_poly[0][0], ee_poly[1][0], ee_poly[2][0], ee_poly[3][0], ee_poly[4][0]]
        y_new = [ee_poly[0][1], ee_poly[1][1], ee_poly[2][1], ee_poly[3][1], ee_poly[4][1]]
        
        number_pixel_x = len(sel_polygon_npy[list_of_bands[0]][0])
        number_pixel_y = len(sel_polygon_npy[list_of_bands[0]])
        
        # min- and max-boundaries are caclulated for every direction
        west = min(x_new) + ((max(x_new) - min(x_new))/number_pixel_x)/2
        east = max(x_new) - ((max(x_new) - min(x_new))/number_pixel_x)/2
        south = min(y_new) + ((max(y_new) - min(y_new))/number_pixel_y)/2
        north = max(y_new) - ((max(y_new) - min(y_new))/number_pixel_y)/2

        # use boundaries to create the coordinates of the xarray 
        ncol_coords = np.linspace(start=south, stop=north, num=len(sel_polygon_npy))
        nrow_coords = np.linspace(start=west, stop=east, num=len(sel_polygon_npy[0]))

        # iterate through the list of bands, extract the selected burn hours and merge them into an xarray dataset
        bands_to_merge = []
        for sel_hour in range_hours:
            band = sel_polygon_npy[list_of_bands[0][sel_hour]] 
            band_to_xr = xarray.DataArray(np.flip(band, axis=0), 
                        coords={'y': ncol_coords,'x': nrow_coords}, 
                        dims=["y", "x"])
            band_to_xr.name = list_of_bands[0][sel_hour]
            band_to_xr = band_to_xr.where(~(band_to_xr==9999),np.nan)
            band_to_xr.rio.write_nodata(np.nan, inplace=True)
            bands_to_merge.append(band_to_xr)
            data_xr = xarray.merge(bands_to_merge)
        
        return data_xr
    
    
    def set_transform_by_hand(self, xarray_ds):
        # Extract the coordinate information
        x_coords = xarray_ds.x.values
        y_coords = xarray_ds.y.values

        # Check the number of pixels in x and y dimensions
        num_pixels_x = len(x_coords)
        num_pixels_y = len(y_coords)

        # Compute pixel width (resolution in x-direction)
        if num_pixels_x > 1:
            pixel_width = (x_coords[-1] - x_coords[0]) / (num_pixels_x - 1)
        else:
            # Define a default pixel width if only one pixel is present
            pixel_width = 0.1  # Adjust this value based on your specific context. Using CRS 4326, 0.1 is the length and height of a pixel

        # Compute pixel height (resolution in y-direction)
        if num_pixels_y > 1:
            pixel_height = (y_coords[-1] - y_coords[0]) / (num_pixels_y - 1)
        else:
            # Define a default pixel height if only one pixel is present
            pixel_height = 0.1  # Adjust this value based on your specific context. Using CRS 4326, 0.1 is the length and height of a pixel

        # Upper left corner coordinates (assuming the dataset is in north-up orientation)
        upper_left_x = x_coords[0] - pixel_width / 2
        upper_left_y = y_coords[0] + pixel_height / 2

        # Construct the geotransform array
        transform = [upper_left_x, pixel_width, 0, upper_left_y, 0, -pixel_height]

        # Use affine transformation for the transform
        affine_transform = Affine(transform[1], transform[2], transform[0],
                                transform[4], transform[5], transform[3])

        to_merge = []
        # Iterate through each DataArray within the Dataset and apply the geotransform
        for var_name in xarray_ds.data_vars:
            da = xarray_ds[var_name]
            # Set the transform manually using write_transform
            da.rio.write_transform(affine_transform, inplace=True)
            to_merge.append(da)
            
        # recreate the merged xarray dataset
        return xarray.merge(to_merge)
    
    def download_era5(self, sel_polygon, start_date, start_hour, end_date, len_fire_sequence,
                      out_directory, comp_names = ['u_component_of_wind_10m','v_component_of_wind_10m','dewpoint_temperature_2m','temperature_2m','surface_pressure', 'total_precipitation']):
        '''
        sel_polygon = geopandas object that contains information about geometry, starting date and ID of a selected burned area 
                     (e.g. one row of a geopandas dataframe with equally sized burned area rectangles)
        start_date = starting date of the investigated fire progression sequence according to the burned area raster
        start_hour = exact hour from which to start the fire progression sequence
        end_date = end date of the investigated fire progression sequence according to the burned area raster
        len_fire_sequence = length of the investigated fire sequence. This corresponds to the hours of the fire propagation interval
        num_sub_images = Defines, in how many sub-images the original image is divided. 
                        Depending on the size of the selected burned area polygon,
                        the resulting sentinel-2 images can be too large to transfer from the google earth engine server to the client side.
                        the original image can therefore be split into several sub-images, which are stitched together again after the download.
                        3 images per row and column = 9 sub-images, 4 images per row and column = 16 sub-images
        crs = defines the crs. Default is crs 4326
        out_path = path to save the sentinel-2 .tif-images
        '''
        # select the google earth engine dataset 
        era5_land = ee.ImageCollection(self.dataset)
        
        # extract the boundary of the selected polygon
        poly_area = list(sel_polygon.geometry.exterior.coords)
        aoi = ee.Geometry.Polygon(poly_area)  
        
        # using the start and the end day of the fire sequence from the input, two variables are created. 
        # The google earth engine image collection is reduced to these specific days.
        start_date_for_filter = start_date
        if start_date != end_date:
            end_date_for_filter = str(pd.to_datetime(end_date)+datetime.timedelta(days=1))[0:10]
        elif np.logical_and(start_date == end_date, int(start_hour) < 12):
            end_date_for_filter = str(pd.to_datetime(end_date)+datetime.timedelta(days=1))[0:10]
        elif np.logical_and(start_date == end_date, int(start_hour) >= 12):
            end_date_for_filter = str(pd.to_datetime(end_date)+datetime.timedelta(days=2))[0:10]
            
        # Using this information, a list is created that conatains all the bands to download     
        range_bands = range(0,((pd.to_datetime(end_date_for_filter) - pd.to_datetime(start_date_for_filter)).days)*24)
        # from this list, the specific hours of the fire sequence are extracted and saved
        range_hours = range(int(start_hour), (int(start_hour)+len_fire_sequence))
        # The google earth engine image collection is reduced to the specified days
        era5_land_sel_poly = (era5_land.filter(ee.Filter.date(start_date_for_filter, end_date_for_filter)))

        # iteratively, all bands of every selected component ("comp_names") are extracted from the image collection, 
        # clipped to the area of interest and saved in a list
        # the band names are also saved in a list
        comp_list = []
        band_names = []
        for idx, comp in enumerate(comp_names):
            
            def sel_component(img):  
                return img.select(comp)
            
            era5_comp = era5_land_sel_poly.map(sel_component)
            
            band_names_comp = [f'{comp}_{m}' for m in range_bands]
            img_of_comp = ee.Image(era5_comp.toBands().rename(band_names_comp))
            img_of_comp = img_of_comp.clipToBoundsAndScale(
                            geometry=aoi, scale=11131.949079327358
                    )

            band_names.append([band_names_comp, idx])
            comp_list.append([img_of_comp, idx])

        # images of all components are transfered from the server (google earth engine) to the client.
        # The different components are downloaded parallelly to increase speed.
        # To be able to assign the results to the correct component of the ERA5 dataset, an index is saved as well
        with concurrent.futures.ThreadPoolExecutor(max_workers=100) as executor:
            single_comp = []
            for era_comp in comp_list:
                single_comp.append(executor.submit(self.get_single_images, ee_image=era_comp[0], index=era_comp[1]))
                
        list_of_results_init = []
        index_list = []
        for future in concurrent.futures.as_completed(single_comp):
            np_arr, index = future.result()
            index_list.append(index)
            list_of_results_init.append(np_arr)       
          
        # Check if the download worked
        if  all(len(x) !=0 for x in list_of_results_init):
            
            # transfrom the individual numpy arrays of every component into an xarray with the same geographical dimensions as the google earth engine image
            # the xarrays and the variable names of every component are saved in a list
            list_comp_arr = []
            list_comp_var = []
            for i in range(len(comp_names)):
                sel_comp_to_xr = list_of_results_init[np.where(np.array(index_list) == i)[0][0]]
                sel_comp_band_names = band_names[i]

                sel_comp_array = self.create_xr_dataset(sel_comp_to_xr, comp_list[i][0], sel_comp_band_names, range_hours)
                list_comp_arr.append(sel_comp_array)
                list_comp_var.append(list(sel_comp_array.data_vars))
            
            # create the directory, in which the ERA5-data of every hour of one fire sequence will be saved
            directory = f"{out_directory}/{i}"
            if not os.path.exists(directory):
                os.makedirs(directory)
            
            # iteratively merge all ERA5-components of one hour of the fire sequence into an xarray dataset
            # and save a raster for every hour of the fire sequence. 
            # In the saved raster image, the exact hour is indicated in the name. 
            for k in range(len(list_comp_var[0])):
                comp_to_merge = []
                for comp in range(len(list_comp_arr)):
                    comp_to_merge.append(list_comp_arr[comp][list_comp_var[comp][k]])
                    
                    meteo_vars = xarray.merge(comp_to_merge)
                    meteo_vars = meteo_vars.rio.write_crs("epsg:4326", inplace=True)
                    
                    # if the number of pixels in y- and/or x-direction, gdal has problems with the transformation
                    # in such cases, the saved raster file won't have any geotransformation attached.
                    # This can be prevented by setting the transformation by hand, which is done here.
                    # Thereto, it is Checked if the number of pixels in x and y exceed 1.
                    # If this isn't the case in either one of the directions, the corresponding function is applied.
                    # Otherwise, nothing is done.
                    if len(meteo_vars.x.values) == 1 or len(meteo_vars.y.values) == 1:
                        meteo_vars = self.set_transform_by_hand(meteo_vars)
                    
                    meteo_vars.rio.to_raster(directory + f"/{sel_polygon['id']}_{start_date}T{str(range_hours[k])}-00-00.tif")
            # if all hours of a fire sequence were saved successfully, this message is printed
            print(f"All ERA5-hours of {sel_polygon['id']} successfully saved")
        else:
            # if it did not work, print the problematic id
            print(f"ERA5-data of {sel_polygon['id']} could not be saved")  