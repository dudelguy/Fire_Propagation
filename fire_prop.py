# import necessary libraries
import pandas as pd
import geopandas as geo
import xarray
import rioxarray as rio
from osgeo import osr, ogr, gdal
import numpy as np
from rasterio import features
from rasterio.mask import mask
from rasterio.enums import MergeAlg
from scipy.interpolate import griddata
from shapely.ops import unary_union
from shapely.geometry import Point, mapping, box
import scipy.ndimage as ndi
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist


class gis_functions:
    
    def __init__(self):
        self = self

    def poly_from_raster(self, input_raster):
        
        '''
        function that creates a polygon from a binary input raster.
        helper function needed to create the .tif rasters of identical size

        input:
        input_raster = raster created during "fire_prog_raster_ident_size". 

        output:
        union_vector_gdf = geodataframe that contains a polygon with the burned area at the date of the input raster
        '''
        
        # get geographical extend (x and y), as well as pixel values of the input raster 
        x, y, time = input_raster.x.values, input_raster.y.values, input_raster.values
        # create meshgrid containing the geographical information of every pixel of the raster
        x, y = np.meshgrid(x, y)
        # assign the geographical information of a pixel to the corresponding pixel value of the raster
        x, y, time = x.flatten(), y.flatten(), time.flatten()

        # create pandas dataframe from this information 
        df = pd.DataFrame.from_dict({'time': time, 'x': x, 'y': y})
        
        # get the position of all burned pixels (= 1) from the raster
        if len(sorted(np.unique(input_raster))) > 1:
            gdf_threshold = sorted(np.unique(input_raster))[1]
        else:
            gdf_threshold = sorted(np.unique(input_raster))[0]
            
        # select all rows of the dataframe that include burned pixels    
        df = df[df['time'] == gdf_threshold]
        
        # create a geodataframe for these pixel. Every pixel is a single geometrical object
        gdf_vector = geo.GeoDataFrame(geometry=geo.GeoSeries.from_xy(df['x'], df['y'], crs=input_raster.rio.crs))
        # draw buffer around every pixel/geometrical object and join the buffered area to create a single polygon. The value 12.5 has to be set to fit the corresponding data
        gdf_vector = gdf_vector.buffer(12.5, cap_style=3)
        total_union = unary_union(gdf_vector)    
        # Create the final GeoDataFrame, containing only one polygon created from the joined pixel
        union_vector_gdf = geo.GeoDataFrame(geometry=geo.GeoSeries(total_union))
        
        # return the final geodataframe 
        return union_vector_gdf       
       

    def create_border_mask(self, binary_mask, radius = 1):
        '''
        function to create a border mask around a given patch of pixel with identical values
        here, all pixels surrounding a burned area of a specific date are selected

        input: 
        binary_mask = binary xarray (here: all burned pixel of a specific date are 1 and non-burned pixel are 0)
        radius = size of the border mask in pixel

        output:
        border_mask = the border mask, indicating all pixels that belong to the specific border
        '''
        # Dilate the binary mask to include neighboring pixels
        radius = radius  # Adjust the radius based on how many pixels you want to include
        dilated_mask = ndi.binary_dilation(binary_mask, iterations=radius)

        # Erode the dilated mask to include only the border values
        border_mask = dilated_mask & ~ndi.binary_erosion(dilated_mask, iterations=1, border_value=1)
        # return the border mask
        return border_mask
       

    def create_raster_coordinates_from_shp(self, input_geometry, res):
        '''
        function used to create the boundaries of an empty raster with the size of the shapefile

        input: 
        input_geometry = geometry for which to create the boundary 
        res = resolution of the raster

        output:
        row_coords = x coordinates of the empty raster
        col_coords = y coordinates of the empty raster
        '''
        # get bounds of the shapefile
        raster_bounds = input_geometry.bounds
        x_min = raster_bounds[0]
        x_max = raster_bounds[2]
        y_min = raster_bounds[1]
        y_max = raster_bounds[3]
        # set resolution of the raster
        raster_res_x = res
        raster_res_y = res

        # calculate the resulting x- and y-coordinates of the raster
        row_coords = np.arange(x_min, x_max, raster_res_x)
        col_coords = np.arange(y_min,y_max, raster_res_y)
        
        return row_coords, col_coords


    def create_raster_from_shp(self, raster_values, input_geometry, res, crs):
        '''
        function that uses the extracted coordinate values of the polygon and
        creates a corresponding raster filled with the specified values

        input:
        raster_values = value to fill the raster with. Can either be a single value to fill the raster with, or an array with values for every pixel of the raster 
        input_geometry = geometry for which to create the boundary 
        res = resolution of the raster
        crs = crs of the raster

        output:
        data_xr = the raster as a rioxarray
        '''
        row_coords, col_coords = self.create_raster_coordinates_from_shp(input_geometry,res)
        
        # create xarray from the interpolated data
        data_xr = xarray.DataArray(raster_values, 
                    coords={'y': col_coords,'x': row_coords}, 
                    dims=["y", "x"])

        # set the crs of the raster to the crs of the input geometry and create a geo xarray
        data_xr = data_xr.rio.write_crs(crs, inplace=True)
        
        return data_xr
       

    def create_new_point(self, fire_id, lat, lon, crs):
        '''
        function to create a geopandas point with specific attributes in the specified crs 

        input: 
        fire_id = string of the fire id, for which the new point is created
        lat = latitude of the point to create
        lon = longitude of the point to create
        crs = crs of the point

        output:
        for_new_point = the new geopandas point
        '''
        # create a new point at the specified coordinates
        geometry = [Point(lon, lat)]
        # create a geopandas dataframe from this point and add the crs, as well as the attributes from the polygon
        for_new_point = geo.GeoDataFrame(geometry=geometry,crs=crs)
        for_new_point["id"] = fire_id
        
        # return the geopandas dataframe that includes this point
        return for_new_point


    def create_and_write_raster(self, output_boundaries, shp_to_write_as_raster, res, crs, out_path):
        '''
        function to create a raster of a burned area, after all pixel of a burned area are assigned to a specific burned date,
        and after the results have been checked for normal fire progression, artifacts, single pixels and islands.
        To ensure that all rasters of all burned areas have the same size, a dataframe with idenitally sized rectangles is needed.
        The resulting raster is saved on disk.

        input:
        output_boundaries = geopandas polygon from the dataset of identically sized rectangles (here, this polygon is geographically linked to the selected burned area) 
        shp_to_write_as_raster = shapefile of the fire progression at the specific firedate
        res = resolution of the output raster
        crs = crs of the output raster. Should have the same crs as the input shp
        out_path = path to save the files 

        output:
        saved .tif file in the specified folder
        '''
        # create empty raster in the size of the output_boundaries (a geopandas polygon)
        template_xr = self.create_raster_from_shp(0, output_boundaries.geometry, res, crs)

        # retrieve coordinates of burned pixels at the specific firedate
        geom_value = list(zip(shp_to_write_as_raster.geometry, [1]))
        # returns an image array with input geometries (in this case: the burned pixels) burned in
        rasterized = features.rasterize(geom_value,
                                        out_shape = template_xr.shape,
                                        transform = template_xr.rio.transform(),
                                        all_touched = True,
                                        fill = 0,   # background value
                                        merge_alg = MergeAlg.add)


        # create xarray from the interpolated data
        data_xr = self.create_raster_from_shp(rasterized, output_boundaries.geometry, res, crs)
        # set all burned pixel to the value 1. All unburned pixel have the value 0
        data_xr = data_xr.where(~(data_xr > 1),1)   
        # write the raster to disk
        data_xr.rio.to_raster(out_path)


class fire_progression(gis_functions):
     
    def __init__(self, fire_det_time, geom):
        
        self.fire_det_time = fire_det_time # name (string) of the column that includes the burn dates (in microseconds). Refers to the geodataframe that contains active fire detection
        self.geom = geom # name (string) of the column that includes the shapefile geometr. Refers to the geodataframe that contains the burned areas
            
    def calc_ignition_point(self, data_xr_clip, ori_cluster_largest):    
        '''
        can be used to approximate the starting pixel of the fire.
        The starting pixel is chosen according to the variance in the distance to the next directly surrounding fire date patches
        The smallest standard deviation indicates smallest variance in the distance, 
        and therefore roughly similar distances in every direction. 
        Accordingly, this point is chosen.

        Input:
        data_xr_clip = xarray with the different burn dates. All artifacts are removed and checked for correct fire progression
        ori_cluster_largest = largest patch of pixels with the oldest burn date. Found with fire_prog_with_iterations

        Output:
        [lat, lon] = List with [Latitude/y, Longitude/x] values of the calculated ignition point
        '''
        # create a binary mask for the largest patch with the oldest fire date
        binary_mask = ori_cluster_largest

        # create border mask for all pixels around the investigated burned area patch
        border_mask = self.create_border_mask(binary_mask, 1)
        
        # extract all pixels surrounding the patch with the oldest fire date
        surrounding_pixels = data_xr_clip.where(border_mask)

        # get values of border pixels (excluding non-fire, i.e. nan, pixels)
        border_without_nan = np.array(~np.isnan(surrounding_pixels))
        # extract coordinates of those border pixels
        border_coordinates = np.argwhere(border_without_nan)
        # extract coordinates of the largest patch with the oldest fire date
        cluster_coordinates = np.argwhere(binary_mask)

        # calculate euclidean distance of a single pixel in the patch with the oldest fire date and every border point
        # this is done for every pixel in the patch
        pairwise_distances = cdist(cluster_coordinates, border_coordinates)
        # create empty list to save the standard deviation of the calcualted distances to the border pixels for every pixel
        std_distances = []

        # iterate through every pixel in the patch and calculate the standard deviations of the distance to the border pixels
        # the values are saved in the empty list
        for i in range(len(pairwise_distances)):
            std_distances.append(np.std(pairwise_distances[i]))

        # get point with smallest standard deviation, i.e. the smalles variance in distance 
        # This ensures that from this pixel, the distance every border pixel is as small as possible
        sel_point_std = cluster_coordinates[np.where(std_distances == np.min(std_distances))[0]]
        
        # extract lon/x and lat/y values of the selected point 
        lon = data_xr_clip.x[sel_point_std[0][1]].values
        lat = data_xr_clip.y[sel_point_std[0][0]].values
        
        # return lon/x and lat/y values of the selected point 
        return [lat, lon]
       

    def fill_poly_with_hulls(self, sel_poly, pot_points_to_sel_fire, crs, res_act_fire = 200):   
        '''
        This function assigns all the pixel of a burned area polygon to their respective burn date. 

        Several steps are performed to fully fill the burned area polygon:
        1- All active fire detection points are buffered according to the resolution of the thermal sensor that was used to detect the active fires 
        2. All active fire detection buffers of an identical date are connected using a convex hull. 
        3. Beginning with the youngest firedate, the convex hulls of every date are "imprinted" on the burned area. 
        This way, pixels that cover multiple burn dates are only assigned to the earliest firedate, since younger, overlapping firedates are overwritten.
        4. The remaining empty pixels are assigned using the k nearest neighbor algorithm. 
        K = 1, and therefore, the pixels are assigned according to the geometrically nearest firedate.
        The nearest neighbour algorithm does not take the shape of the polygon into account, which is why artifacts can occur. 
        They are removed using the functions: rem_artifacts and check_for_islands

        The input are two geopanda-datasets, the output is a rioxarray-raster, in which every pixel is either assigned to a specific burn date or to no-fire

        input:
        sel_poly = the selected burned area polygon (extracted from a geopandas dataframe)
        pot_points_to_sel_fire = geopandas dataframe with all active fire detections assigned to the investigated burned area polygon. 
                                If VIIRS data is used, the output from "find_viirs_for_effis" can be used
        res_act_fire = resolution of the thermal sensor used for active fire detection
        crs = the crs of the ouput

        ouput:
        data_xr_clip = a rioxarray raster with fire/no-fire pixels. Fire piels are labelled according to their firedate (in microseconds)
        '''
        # add new column to active fire dataset, in which each active fire detection point is buffered according to the resolution of the thermal sensor
        pot_points_to_sel_fire["buffered_points"] = pot_points_to_sel_fire.buffer(res_act_fire)
        # save all unique fire dates in a variable (used to iterate through every individual burn date)
        unique_timestamps = np.unique(pot_points_to_sel_fire[self.fire_det_time])

        # empty list to save the convex hulls
        hulls = []
        # iterate through every distinct burn date. 
        # Select all points of one burn date and construct the convex hull for this specific date.
        # Each hull is saved in the empty list, together with the firedate
        for i in range(len(unique_timestamps)):
            sel_points = pot_points_to_sel_fire[pot_points_to_sel_fire[self.fire_det_time]==unique_timestamps[i]]
            hull_of_points = unary_union(sel_points.buffered_points.tolist()).convex_hull
            hulls.append([unique_timestamps[i], hull_of_points])
            
        # construct a geopandas dataframe from the convex hulls. 
        # each row contains a hull for a specific burn date
        df_hulls = pd.DataFrame.from_records(data=hulls, columns=['timestamp','geometry'])
        gdf_hulls = geo.GeoDataFrame(data=df_hulls, geometry=df_hulls['geometry'], crs=crs)

        xr_to_fill = self.create_raster_from_shp(np.nan, sel_poly.buffer, 25, crs=crs)
        xr_to_fill = xr_to_fill.rio.clip(geometries=[sel_poly[self.geom]], crs=crs)

        nrow_coords = xr_to_fill.coords['x'].values
        ncol_coords = xr_to_fill.coords['y'].values

        # iterate through the firedates, beginning with the youngest firedate
        for i in reversed(range(len(unique_timestamps))):
            
            # create xarray from the interpolated data.
            # this array is completely filled with the selected firedate
            xr_sel_timestamp = xarray.DataArray(unique_timestamps[i], 
                        coords={'y': ncol_coords,'x': nrow_coords}, 
                        dims=["y", "x"])
            # set the crs and transform it to a rio xarray
            xr_sel_timestamp = xr_sel_timestamp.rio.write_crs(crs, inplace=True)
            # clip the raster to fit the convex hull of the specific burn date
            xr_clip_sel_timestamp = xr_sel_timestamp.rio.clip(geometries=[gdf_hulls.iloc[i]["geometry"]], crs=crs, drop = False)
            # clip the convex hull to the burned area 
            xr_clip_sel_timestamp = xr_clip_sel_timestamp.rio.clip(geometries=[sel_poly[self.geom]], crs=crs)
            
            # select all pixels that, after clipping to the burned area, still have a firedate
            to_fill = np.argwhere(~np.isnan(np.array(xr_clip_sel_timestamp)))

            # iterate over every pixel and fill all pixels of the initially empty raster with the corresponding firedates
            for filling in range(len(to_fill)):
                xr_to_fill[to_fill[filling][0],to_fill[filling][1]] = xr_clip_sel_timestamp[to_fill[filling][0],to_fill[filling][1]]

        # create meshgrid of the y and x values (lat/lon)
        x_grid, y_grid = np.meshgrid(nrow_coords, ncol_coords)
        
        # extract and store pixel values and pixel coordinates of all pixel that are filles with a firedate
        point_values = np.concatenate(xr_to_fill.values)
        point_coords = list(zip(np.concatenate(x_grid),np.concatenate(y_grid)))
        idx_nan = np.where(np.isnan(point_values))[0]
        idx_0 = np.where(point_values == 0)[0]
        idx_to_rem = np.concatenate((idx_nan.astype(int),idx_0.astype(int)))
        idx_to_keep = np.arange(0, len(point_values), 1)
        idx_to_keep = np.delete(idx_to_keep, idx_to_rem)
        point_values = np.delete(point_values, idx_to_rem)
        point_coords_keep = []
        for i in range(len(idx_to_keep)):
            point_coords_keep.append(point_coords[idx_to_keep[i]])  

        # create an empty grid, in which the formerly extracted pixel are assigned to a value 
        # all remaining, empty pixels are filled with a firedate by applying the nearest neighbor approach
        # i.e. they are filled with the firedate of the geometrically nearest pixel that was initially assigned to a firedate
        grid_date = griddata(point_coords_keep,point_values,(x_grid, y_grid), method="nearest")

        # create xarray from the interpolated data
        data_xr = xarray.DataArray(grid_date, 
                    coords={'y': ncol_coords,'x': nrow_coords}, 
                    dims=["y", "x"])

        # set the crs and make it rioxarray
        data_xr = data_xr.rio.write_crs(crs, inplace=True)
        # clip the final result to the burned area polygon
        data_xr_clip = data_xr.rio.clip([sel_poly[self.geom]], crs=crs)

        return data_xr_clip


    def fire_prog_raster_ident_size(self, filled_raster, selected_rectangle, res, crs, dates_str, out_path):   
        '''
        saves binary rasters of identical size to disk. 
        One file is created for every fire date of the chosen, iteratively increasing fire size with every new fire date.
        All files of one fire have the same height and width 

        input:
        filled_raster = raster, in which the pixels are filled with zeroes (no fire) or distinct dates of fire propagation (here, dates are given in microseconds)
        selected_poly = used to link the saved file to a fire ID (here, the fire ID is the same as in the EFFIS database)
        output_path = path to store the saved files

        output:
        .tif file saved in the folder of the selected output path
        '''

        # get all distinct dates of fire progression
        # these dates are based on active fire detections converted to microseconds
        unique_dates = np.unique(filled_raster)[~np.isnan(np.unique(filled_raster))]

        # create an empty raster with the same size as the raster that is filled with the dates
        raster_fire_prog = filled_raster.copy()
        raster_fire_prog = raster_fire_prog.where(np.isnan(raster_fire_prog),0)
        raster_fire_prog = raster_fire_prog.where(~np.isnan(raster_fire_prog),0)

        # iterate through the distinct dates and fill the empty raster accordingly -> 0 = no fire, 1 fire at date i
        # For every fire date, a file is saved to disk
        for i in range(len(unique_dates)):
            # create a binary mask to identify all pixels of a distinct date
            binary_mask = (filled_raster == unique_dates[i])  
            
            # use binary mask to create a new raster with 0 (no fire at date i) and 1 (fire at date i)
            fire_prog = filled_raster.where(binary_mask)    
            fire_prog = fire_prog.where(np.isnan(fire_prog),1)    
            fire_prog = fire_prog.where(~np.isnan(fire_prog),0)    
            
            # newly burned areas are added during the loop. 
            # the resulting raster is getting consecutiely filled with fire pixels, keeping older burned areas as well
            raster_fire_prog = np.add(raster_fire_prog, fire_prog)
            
            # create polygon of this raster. Needed for create_and_write_raster
            fire_prog_shp = self.poly_from_raster(raster_fire_prog)
            # use function to create and save the raster
            out_path_iter =f'{out_path}/{selected_rectangle["id"]}_{i}_{dates_str[i]}.tif'
            self.create_and_write_raster(selected_rectangle, fire_prog_shp, res, crs, out_path_iter)
            

class active_fire_detection():
    
    def __init__(self, fire_detection_time, fire_detection_date, lat, lon, fire_id, fire_start_db, fire_end_db):
        '''
        All of these variables are strings with the names of the respective columns in the active fire detection .csv file
        '''
        self.acq_time = fire_detection_time
        self.acq_date = fire_detection_date     
        self.lat = lat
        self.lon = lon
        '''
        All of these variables are strings with the names of the respective columns in the burned area geopandas dataframe
        '''
        self.id = fire_id 
        self.start_fire = fire_start_db
        self.end_fire = fire_end_db   


    def select_subset_using_id(self, ori_dataset, subset_of_ids):
        '''
        function to select a subset of burned area polygons from a dataset based on a list of ids

        input: 
        ori_dataset = the dataset of burned area polygons. Needs to have a column with unique ids for each polygon
        subset_of_ids = the subset of ids to select from the complete dataset

        output:
        updated_dataset = the selected subset
        '''
        # extract indices at which the ids of the subset can be found in the complete dataset
        idx_selection = ori_dataset[self.id].isin(subset_of_ids)
        idx_true = np.where(idx_selection==True)[0]
        # extract the subset based on these indices
        updated_dataset = ori_dataset.iloc[idx_true]
        
        return updated_dataset

    def index_containing_substring(self, list_with_act_fire_files, substring):
        '''
        the active fire detections are downloaded in a country-specific .csv file (from e.g. https://firms.modaps.eosdis.nasa.gov/download/create.php)
        this function is used to find these country specific VIIRS files in a folder that contains all country-files

        input:
        list_with_act_fire_files = list that contains the names of .csv-files of all countries
        substring = abbreviation of the investigated country

        output: 
        i = name of the country specific .csv file 
        '''
        # iterate throuhg the list of files and return the wanted filename as soon as it is found
        for i, s in enumerate(list_with_act_fire_files):
            if substring in s:
                return i
    
    def modify_viirs(self, sel_country_and_year, crs="EPSG:4326"):
        '''
        function that transfroms a dataframe with the active fire detections into a geopandas dataframe
        Also adds a column "datetime", in which the aquisition time of the active fire detection is combined with the acquisition date 

        input: 
        sel_country_and_year = geodataframe with active fire detections for a specific country and a specific year

        output: 
        active_fire_gdf = an updated geopandas dataframe 
        '''
        lon = self.lon
        lat = self.lat
        
        # select the column that contains the acquisition time of the active fire detection
        time1 = sel_country_and_year[self.acq_time].astype(str)
        # iterate over every entry and add 0s to unify the time strings 
        for t in range(len(time1)):
            lenStr = len(time1[t])
            if lenStr == 4:
                continue
            elif lenStr == 1:
                time1[t] = '000' + time1[t]
            elif lenStr == 2:
                time1[t] = '00' + time1[t]
            elif lenStr == 3:
                time1[t] = '0' + time1[t]

        # combine acquisition date and acquisition time transform it into a pandas date and add it as a new column to the dataset 
        sel_country_and_year["datetime"] = pd.to_datetime(sel_country_and_year[self.acq_date].astype(str) + ' ' + time1)

        # tranform the dataframe into a geopandas dataframe
        active_fire_gdf = geo.GeoDataFrame(
            sel_country_and_year, geometry=geo.points_from_xy(sel_country_and_year[lon], sel_country_and_year[lat]), crs=crs
        )
        # return this geopandas dataframe
        return active_fire_gdf
    
    

    def find_viirs_for_effis(self, selected_poly, active_fire_data, day_thresh):
        
        '''
        find all active fire detections that are linked to the selected burned area polygon

        inpout: 
        selected_poly = the selected burned area polygon, already buffered according to the wanted bufferzone
        active_fire_data = the ouput of modify_viirs 
                    -> a geopandas dataframe with active fire detections from the year and the country of the selected burned area polygon
        day_thresh = oftentimes, the date of the burn for the burned area polygon does not cover the complete duration of the fire. 
                    this parameter can be used to include active fire detections from +- x days 

        output: 
        sel_fire_points = a geopandas dataframe that includes all active fire detections for the selected burned area
        '''
        
        # sets search radius for the date of active fire detections
        # uses columnds "LASTUPDATE" and "FIREDATE" from the EFFIS dataset and adds/substracts the additionaly number of days
        # the investigated dataset of burned areas might have different names for the first and last fire dates. 
        # -> change "LASTUPDATE" and "FIREDATE" accordingly
        top_date = pd.to_datetime(selected_poly[self.start_fire]) + pd.Timedelta(days=day_thresh)
        low_date = pd.to_datetime(selected_poly[self.end_fire]) - pd.Timedelta(days=day_thresh)

        # extract all active fire detections inside of the buffered zone of the selected burned area polygon
        sel_fire_points = active_fire_data[active_fire_data.within(selected_poly["buffer"])].copy()

        # select all active fire detections that are inside of the buffer zone AND the inside of the selected time span
        fire_idx = np.where(np.logical_and(sel_fire_points["datetime"]<=top_date, sel_fire_points["datetime"]>=low_date))[0]
        sel_fire_points = sel_fire_points.iloc[fire_idx]

        # return these active fire detections as a geopandas dataframe
        return sel_fire_points
    
        
    def date_to_microseconds(self, sel_fire_points_updated): 
        '''
        for calculations, dates are converted to seconds elapsed since reference time (midnight 1970 UTC)

        input: 
        sel_fire_points_updated = dataframe, which contains the active fire detections that are connected to the chosen fire polygon

        output: 
        sel_fire_points_updated = updated geopandas dataframe with active fire observations 
                                -> a new column "time_in_seconds" is added to the dataframe, representing the seconds elapsed since midnight 1970 UTC

        '''
        # from the datetime, calculate the seconds elapsed since the reference time (midnight 1970 UTC)
        microseconds = []
        for j in range(len(sel_fire_points_updated)):
            microseconds.append(sel_fire_points_updated.iloc[j]["datetime_grouped"].timestamp()) 

        # add the seconds to the original dataframe
        # this is used to automatically select all viirs point of a fire and remove wrong dates   
        sel_fire_points_updated["time_in_seconds"] = microseconds    
        
        return sel_fire_points_updated


    def combine_sat_paths(self, sel_fire_points, dif_hours, dif_minutes):   
        '''
        At European latitudes, the VIIRS instrument has a repetition rate of ~12 hours. Nevertheless, successive satellite paths 
        sometimes observe the same active fire with an offset of 1-2 hours. This function is used to assign a single datetime to these observations

        Input:
        sel_fire_points: output of the find_viirs_for_effis function 
                                -> a geopandas dataframe that includes all active fire detections associated with a selected burned area  
        dif_hours: defined offset in hours
        dif_minutes: defined offset in minutes

        Output:
        pot_points_to_sel_fire_updated: the sel_fire_points geopandas dataframe with an new datetime column (datetime_grouped). 
                                        The new column assigns only a single datetime to observations from successive satellite paths
        '''
            
        # create copy to update
        sel_fire_points_updated = sel_fire_points.copy()
        
        # Calculate time differences between consecutive rows
        time_diff = sel_fire_points_updated['datetime'].diff()
        # Define the time threshold using dif_hours and dif_minutes (we used 1 hour and 41 minutes)
        time_threshold = pd.Timedelta(hours=dif_hours, minutes=dif_minutes)
        # Filter for rows with time difference less than or equal to the threshold
        rows_with_larger_dif = np.where(time_diff > time_threshold)[0]    
        
        # check if any successive satellite path observations exist
        if len(rows_with_larger_dif)>0:   
            # if observations from successive satellite paths exist, create a new array that can be added to the dataframe
            grouped_datetime = np.array([max(sel_fire_points_updated['datetime'])] * len(sel_fire_points_updated))
            
            # iterate over all indices that indicate time differences greater than the threshold
            # all indices in below the selected index are filled with the youngest datetime of this date range
            rows_to_group = 0
            for w in range(len(rows_with_larger_dif)):
                grouped_datetime[rows_to_group:rows_with_larger_dif[w]] = [sel_fire_points_updated['datetime'].iloc[rows_to_group]]
                rows_to_group = rows_with_larger_dif[w]
            
            # the new array "datetime_grouped" is added to the geopandas dataframe    
            sel_fire_points_updated["datetime_grouped"] = grouped_datetime
        else:
            # if there are no observations from successive satellite paths, a new column "datetime_grouped" is still added to the geopandas dataframe  
            # it has the same values as the column "datetime"
            sel_fire_points_updated["datetime_grouped"] = sel_fire_points_updated.iloc[0]['datetime']
        
        
        sel_fire_points_updated = self.date_to_microseconds(sel_fire_points_updated)    
    
        # the updated geopandas dataframe is returned
        return sel_fire_points_updated    
        

class artifacts(gis_functions):
    
    def __init__(self, fire_det_time):
        self.fire_det_time = fire_det_time # name (string) of the column for the burn dates (in microseconds), included in the geodataframe of active fire detection 
        
    def remove_single_pixel(self, fire_prog_array):
        '''
        function to integrate single pixels of a specific fire date into their surrounding.
        this follows the assumption that single pixel are artifacts of the algorithm used to fill the polygon. 

        input: 
        fire_prog_array = xarray after filling it with the different fire propagation dates. I.e., every pixel that belongs to a burned area is filled with its specific burn date (in seconds).
                          In this xarray, individual pixel can have isolated burn dates that do not match that of any neighbour. These are integrated their surounding

        output:
        fire_prog_updated = updated xarray, in which the isolated single pixel are now either assigned to their next closest burn date or, if no burned pixel is in a rook-like neighbourhood, as NA
        '''

        # the output xarray is created from the original xarray
        fire_prog_updated = fire_prog_array.copy()    
        # all unique burn dates are extracted (burn dates are in seconds)
        unique_dates = np.unique(fire_prog_updated)[~np.isnan(np.unique((fire_prog_updated)))]   
        
        # iteration over every unique burn date 
        for dates in unique_dates:
            
            # change every pixel to 0, if it is different from the selected date 
            data_transformed = fire_prog_updated.where(fire_prog_updated == dates, 0)
            # scipys label function is used to identify all consecutive areas of the selected date (rook neighbourhood)
            labeled_array, num_features = ndi.label(data_transformed)

            # a second loop iterates over every area identified with the ndi.label function
            # since ndi.label starts its labelling with 1, the loop aso starts with 1
            for i in range(1,(num_features+1)):
            
                # create a mask for the selected area     
                binary_mask = (labeled_array == i)
        
                # check if the selected mask only contains a single pixel
                if len(np.argwhere(binary_mask)) == 1:
                    surrounding_values = []
                    
                    # if the selected mask really contains an isolated pixel, it is selected and its x and y positions in the xarray are saved 
                    single_pixel_labels = np.argwhere(binary_mask)[0]
        
                    x = single_pixel_labels[0]
                    y = single_pixel_labels[1]
            
                    # All values of the surrounding pixels are selected
                    # If clauses are necessary, since the selected pixel might lie at the border of the raster
                    if (x-1) >= 0:
                        surrounding_values.append(fire_prog_updated[x-1,y].values.item())
                    if (y-1) >= 0:
                        surrounding_values.append(fire_prog_updated[x,y-1].values.item())
                    if (y+1) < len(labeled_array[0]):
                        surrounding_values.append(fire_prog_updated[x,y+1].values.item())
                    if (x+1) < len(labeled_array):
                        surrounding_values.append(fire_prog_updated[x+1,y].values.item())
    
                    # check if any surrounding pixel is filled with a burn date. 
                    # if yes, it is filled with the next highest date
                    # if not, it is filled with na
                    if len(np.array(surrounding_values)[~np.isnan(np.array(surrounding_values))]) > 0:
                        fire_prog_updated[x,y] = min(np.array(surrounding_values)[~np.isnan(np.array(surrounding_values))])
                    else:
                        fire_prog_updated[x,y] = np.nan    

        return fire_prog_updated
    

    def check_for_islands(self, fire_prog_array):
        '''
        function to check for unregularities inside of burned areas.
        more specifically, if a burned area of a specific date is inside the burned area of another date

        input:
        fire_prog_array = xarray after filling it with the different fire propagation dates. I.e., every pixel that belongs to a burned area is filled with its specific burn date (in seconds).
                          In this xarray, unregularities can occur, concerning burn dates that do not match the burn dates of their neighbours. These "islands" are integrated into their surrounding.

        output:
        fire_prog_updated = new xarray, in which the affected areas are removed 
        '''
        
        # the output xarray is created from the original xarray
        fire_prog_updated = fire_prog_array.copy()    
        # all unique burn dates are extracted (burn dates are in seconds)
        unique_dates = np.unique(fire_prog_updated)[~np.isnan(np.unique((fire_prog_updated)))]   
        
        # check for the number of individual burn dates
        # if there is are only two burn date, burned "islands" cannot be present and this check is obsolete
        if len(unique_dates)>2:
            
            # iterate over the individual burn dates
            # we always assume that the first date was filled correctly -> start iteration with the second fire date
            for dates in range(1,(len(unique_dates)-1)):
                
                # select the date to be checked
                inv_date = unique_dates[dates]
                # set all pixels from different dates to zero                
                data_transformed = fire_prog_updated.where(fire_prog_updated == inv_date, 0)
                # select all younger burn dates
                younger_dates = unique_dates[(dates+1):(len(unique_dates)+1)]
                # select all older burn dates
                older_dates = unique_dates[0:(dates)]
                
                # scipys label function is used to identify all consecutive areas of the selected date (rook neighbourhood)
                labeled_array, num_features = ndi.label(data_transformed)

                # iterate over every area of the selected burn date (identified with the ndi.label function)
                # since ndi.label starts its labelling with 1, the loop aso starts with 1
                for i in range(1,(num_features+1)):
                    
                    # create a mask for the selected area     
                    binary_mask = (labeled_array == i)

                    #create the border mask for the selected date
                    border_mask = self.create_border_mask(binary_mask)

                    # Use this mask to extract all surrounding pixels
                    surrounding_pixels = np.array(fire_prog_updated)[border_mask]
                                
                    # get minimum date from surrounding pixels. Na values are included in this array
                    unique_dates_surrounding_withNA = np.unique(surrounding_pixels)
                    # get minimum date from surrounding pixels. Na values are excluded in this array
                    unique_dates_surrounding = np.unique(surrounding_pixels)[~np.isnan(np.unique((surrounding_pixels)))]
                    
                    # check, if the selected area is an island inside of another burned area (no NAN values in border mask), 
                    if len(unique_dates_surrounding) == 1 and len(unique_dates_surrounding_withNA) == 1:
                        # we started with the second burn date and we used convex hulls to create consecutive burn patches
                        # therefore, if a patch of pixels of the second burn date is inside another burn date, this other burn date has to be younger.
                        # the resulting burn date has therefore no direct connection to an older burn date. 
                        # This is unrealisticfire propagation, and we remove these areas by replacing the affected pixels with the next younder burn date
                        fire_prog_updated = fire_prog_updated.where(labeled_array != i,min(unique_dates_surrounding))    
                    
                    # if it is at the border of the final burned area (NAN values in border mask),
                    elif len(unique_dates_surrounding) == 1 and len(unique_dates_surrounding_withNA) == 2:
                        # if the burned patch is at the border of the burned area and has connection to an older burned area, nothing is changed
                        if np.any(np.isin(unique_dates_surrounding,older_dates)):
                            fire_prog_updated = fire_prog_updated
                        
                        # if the burned patch is at the border of the burned area and has no connection to an older burned area, 
                        # all pixels of this patch are replaced by the youngest burn date surrounding the patch
                        else:
                            fire_prog_updated = fire_prog_updated.where(labeled_array != i,max(unique_dates_surrounding)) 
                            
                    # if it is completely isolated from other burned areas,
                    elif len(unique_dates_surrounding) == 0 and len(unique_dates_surrounding_withNA) == 1:
                        # nothing is changed
                        fire_prog_updated = fire_prog_updated
                        
                    # or if the surrounding dates all belong to newer burn dates
                    elif np.all(np.isin(unique_dates_surrounding,younger_dates)):    
                        # if this is true after all the former checks, all pixels of the selected patch are replaced by the next younger burn date
                        fire_prog_updated = fire_prog_updated.where(labeled_array != i,min(unique_dates_surrounding))                
            
            # return the updated arrays 
            return fire_prog_updated   
        else:
            # nothing needs to be changed -> return the original array
            return fire_prog_array
        
        
    def rem_artifacts(self, fire_prog_array, pot_points_to_sel_fire, distance_thresh):   
        '''
        To link a pixel to a firedate, the nearest neighbor algorithm is used on the rectangular xarray raster.
        This does not consider the actual shape of the polygon.
        As a result, in specific cases and for specific burned area shapes, 
        patches of the polyon are assigned without any connection to an active fire detection measurement

        input:
        fire_prog_array = xarray after filling it with the different fire propagation dates. I.e., every pixel that belongs to a burned area is filled with its specific burn date (in seconds).
                          In this xarray, irregularities at the borders of the final burned area can occur.

        pot_points_to_sel_fire = active fire detection points that are assigned to the burned area

        distance_thresh = distance in meters. If no actove fire point is in the selected distance of the burned area patch, the values of this patch are replaced  


        output:
        fire_prog_updated = new xarray, in which the affected areas are correctly assigned 
        '''
        
        # the output xarray is created from the original xarray
        fire_prog_updated = fire_prog_array.copy()    
        # all unique burn dates are extracted (burn dates are in seconds)
        unique_dates = np.unique(fire_prog_updated)[~np.isnan(np.unique((fire_prog_updated)))]   

        # select all fire dates of the VIIRS detections (in seconds)
        point_values = pot_points_to_sel_fire[self.fire_det_time].values
        # extract the coordinates of all VIIRS point measuements 
        point_coords = list(zip(pot_points_to_sel_fire.geometry.x, pot_points_to_sel_fire.geometry.y))

        # iterate over all individual burn dates
        for dates in unique_dates:
            
            # select coordinates of th selected date    
            point_coords_inv_date = np.array(point_coords)[(np.where(point_values==dates)[0])]        
            # set all pixels from different dates to zero                
            data_transformed = fire_prog_updated.where(fire_prog_updated == dates, 0)
            # scipys label function is used to identify all consecutive areas of the selected date (rook neighbourhood)
            labeled_array, num_features = ndi.label(data_transformed)            

            # iterate over every area identified with the ndi.label function
            # since ndi.label starts its labelling with 1, the loop aso starts with 1
            for i in range(1,(num_features+1)):
                    
                # create a mask for the selected area         
                binary_mask = (labeled_array == i)
                
                # select coordinates of pixels belonging to the selected area
                cluster_coords_xy = np.argwhere(binary_mask)

                # create list that saves coordinate pairs for all pixels of the selected area
                cluster_coords_latlon = []
                for inv_cluster_pixel in range(len(cluster_coords_xy)):
                    lon = fire_prog_array[cluster_coords_xy[inv_cluster_pixel][0],cluster_coords_xy[inv_cluster_pixel][1]].x.values.item()
                    lat = fire_prog_array[cluster_coords_xy[inv_cluster_pixel][0],cluster_coords_xy[inv_cluster_pixel][1]].y.values.item()
                    
                    cluster_coords_latlon.append((lon, lat))
                
                # scipys distance matrix calculates distances for all points of the selected area and the VIIRS detections
                dist_to_points = distance_matrix(cluster_coords_latlon, point_coords_inv_date)
                
                # check if distance of all points to any VIIRS detection is greater than the the distance_thresh
                if np.all(dist_to_points.flatten() > distance_thresh):
                    # if this is the case, create the border mask for the selected date
                    border_values = self.create_border_mask(binary_mask)
                    # get unique values of all non-na border pixels
                    surrounding_pixels = np.array(fire_prog_updated)[border_values]
                    surrounding_pixels = np.unique(surrounding_pixels)[~np.isnan(np.unique(surrounding_pixels))]
                    
                    # if there is at least one non-na boder pixel, the pixels of the selected area are replaced by it.
                    # if there multiplte dates at the border of the selected area, the next youngest date is used 
                    if len(surrounding_pixels) > 0:
                        fire_prog_updated = fire_prog_updated.where(~binary_mask, min(surrounding_pixels))
                    else:
                        fire_prog_updated = fire_prog_updated.where(~binary_mask, np.nan)  
                    
                else:
                    continue
                
        return fire_prog_updated    


class check_fire_progression(gis_functions):
    
    def __init__(self):
        self 

    def new_minima(self, fire_prog_array, labeled_array, investigated_feature):
        '''
        Helper function to check for "normal" fire progression 
        Neighbors of selected burn date are checked for continous fire progression.
        Afterwards, pixel of the selected burn date are overwritten with the next younger burn date of its direct neighborhood

        input:
        fire_prog_array = xarray, in which every pixel is assigned to a different burn date (or no burn). All artifacts are already removed
        labeled_array = first output of the ndi.label() function
        investigated_feature = specific number of the investigated features returned by the ndi.label() function

        output: 
        clipped_dates_new_min = new xarray, in which pixels of the investigated burn date are overwritten with the next younger burn date

        Also returns one of two number:
        1) if fire date of the selected burned area has the next younger fire date as a neighbor, 1 is returned 
        2) if fire date of the selected burned area has an older date as neighbor, 0 is returned 

        results are used to monitor stepwise, if the fire progresses in a natural way (from older dates to newer ones)
        '''
        # get dates of pixels of burned areas
        old_min_date = np.unique(fire_prog_array)
        # create binary mask for the specific burned area patch
        binary_mask = (labeled_array == investigated_feature)
        # create border mask for all pixels around the investigated burned area patch
        border_mask = self.create_border_mask(binary_mask, radius = 1)
        # Use border mask to extract all surrounding pixels
        surrounding_pixels = fire_prog_array.where(border_mask)
                    
        # get minimum date from surrounding pixel
        min_date_surrounding_pixels = min(np.unique(surrounding_pixels))
  
        # check, if the surrounding pixels include the next younger date          
        if min_date_surrounding_pixels == old_min_date[1]:
            # set dates of the selected burned area patch to the youngest date found in the border pixel
            clipped_dates_new_min = fire_prog_array.where(labeled_array != investigated_feature,min_date_surrounding_pixels) 
            # if the new date is also the next younger date, return the value 1    
            return clipped_dates_new_min, 1
        else:
            # if the new date is not the next younger date of the investigated burned area, 
            # return the value 0 (no normal fire progression observed)       
            return fire_prog_array, 0 
        

    def fire_progression(self, fire_prog_array):
        '''
        Helper function to check for "normal" fire progression.
        Uses new_minima function on all areas from the selected fire date and return all 
        fire date patches that show normal fire progression (i.e. next younger fire date in direct neighborhood)

        input:
        fire_prog_array = xarray with the different burn dates. All artifacts are already removed

        output:
        clipped_dates_prog = updated fire_prog_array, in which the oldest fire date is replaced by the next younger date. 
                             Can be used iteratively, to check for correct fire progression over multiple iterations
        cluster_mask_with_correct_progression = list of all patches that show normal fire progression for the investigated fire date
        '''
        # copy of the original array, which will be updated and returned
        clipped_dates_prog = fire_prog_array.copy()

        # initialize empty list to store the mask of the correct progression
        cluster_mask_with_correct_progression = []
        # set all pixels except those of the earliest fire date to 0               
        data_transformed = clipped_dates_prog.where(clipped_dates_prog == min(np.unique(clipped_dates_prog)), 0)
        # label all the individual patches that belong to the specific fire date
        labeled_array, num_features = ndi.label(data_transformed)
            
        # iterate over every patch identified with the ndi.label function
        # since ndi.label starts its labelling with 1, the loop aso starts with 1
        for i in range(1,(num_features+1)):
                
            # apply new_minima function. If is_progressing is 1, the fire is progressing correctly
            clipped_dates_prog, is_progressing = self.new_minima(clipped_dates_prog,labeled_array,i)
            if is_progressing == 1:
                # if the fire progression worked, an xarray is saved, 
                # in which all pixel of the investigated patch are True, 
                # while all remaining pixel are False
                cluster_mask_with_correct_progression.append(labeled_array == i)
            
        return clipped_dates_prog, cluster_mask_with_correct_progression     


    def fire_progression_v2(self, fire_prog_array):
        '''
        Helper function to check for "normal" fire progression.
        Uses new_minima function on all areas from the selected fire date and return all 
        fire date patches that show normal fire progression (i.e. next younger fire date in direct neighborhood)

        input:
        fire_prog_array = xarray with the different burn dates. All artifacts are already removed

        output:
        clipped_dates_prog = updated fire_prog_array, in which the oldest fire date is replaced by the next younger date. 
                                Can be used iteratively, to check for correct fire progression over multiple iterations
        cluster_mask_with_correct_progression = list of all patches that show normal fire progression for the investigated fire date
        '''
        # copy of the original array, which will be updated and returned
        clipped_dates_prog = fire_prog_array.copy()

        # initialize empty list to store the mask of the correct progression
        cluster_mask_with_correct_progression = []
        # get all burn dates 
        burn_dates = np.unique(clipped_dates_prog)
        # set all pixels except those of the earliest fire date to 0               
        data_transformed = clipped_dates_prog.where(clipped_dates_prog == min(burn_dates), 0)
        # label all the individual patches that belong to the specific fire date
        labeled_array, num_features = ndi.label(data_transformed)
        # set all pixels except those of the second earliest fire date to 0               
        burn_mask_2 = clipped_dates_prog.where(clipped_dates_prog == burn_dates[1], 0)
        # label all the individual patches that belong to the specific fire date
        labeled_array_next_date, num_features_next_date = ndi.label(burn_mask_2)
        
        # iterate over every patch identified with the ndi.label function
        # since ndi.label starts its labelling with 1, the loop aso starts with 1
        
        fire_prop_cluster = []
        for i in range(1,(num_features+1)):
                
            feature_mask = labeled_array == i   
            clump_border = self.create_border_mask(feature_mask) 
            
            for k in range(1, (num_features_next_date+1)): 
                    
                binary_mask = (labeled_array_next_date == k)    
                if np.any((binary_mask) & clump_border) == True:
                    
                    clipped_dates_prog = clipped_dates_prog.where(labeled_array != i,burn_dates[1])
                    cluster_mask_with_correct_progression.append(labeled_array == i)
                    fire_prop_cluster.append(k)
        
        fire_prop_cluster = np.unique(fire_prop_cluster)
        if len(fire_prop_cluster) == num_features_next_date:    
            return clipped_dates_prog, cluster_mask_with_correct_progression     
        else: 
            return clipped_dates_prog, []

    def fire_prog_iterations(self, first_cluster, fire_prog, cluster_with_correct_progression, number_of_iterations):
        '''
        helper function to check for "normal" fire progression
        executes the loop that checks for normal fire progression over n iterations

        all input values are directly given in the fire_prog_with_iterations function
        '''
        # loop over the remaining fire dates and check after every iteration if fire progress is normal.
        # If it is not, leave the loop
        iter = 0
        while iter < (number_of_iterations):
            if len(cluster_with_correct_progression)>0:
                fire_prog, cluster_with_correct_progression = self.fire_progression_v2(fire_prog)
            else:
                break
            iter = iter + 1
        
        # initialize empty list to save future pixel patches that progress normally
        ori_cluster_list = []
        
        if len(cluster_with_correct_progression)>0:
            merged_array = np.any(cluster_with_correct_progression, axis=0)
            # check if the results overlap with the initial check for fire progression
            # if it does, save the corresponding starting patch
            for ori_clust in first_cluster:  
                ori_cluster_list.append(np.logical_and(ori_clust == True, merged_array == True))       
            # if the results do not overlap with the inital starting dates, return an empty list    
            if ~ori_cluster_list[0].any(): 
                ori_cluster_list = []
                
        return fire_prog, ori_cluster_list


    def fire_prog_with_iterations(self, fire_prog_array, number_of_max_iterations):
        '''
        Final function to check for "normal" fire progression.
        Uses the helper function fire_progression iteratively to monitor the fire progression over a specified number of fire dates

        input:
        fire_prog_array = xarray with the different burn dates. All artifacts are already removed
        number_of_max_iterations = number of consecutive dates for which fire propagation is checked

        output:
        fire_prog = xarray with burned area dates after the specified number_of_max_iterations. 
        first_cluster/ori_cluster_list = if normal fire progression was observed, 
                                         this list includes a binary mask, in which the all pixel of the oldest burn dates are True.
                                         if no normal fire progression was observed, 
                                         this list is empty
        '''
        # number of unique fire dates, excluding non-burned pixels (i.e. na values)
        data_differences_init_fire = len(np.unique((fire_prog_array))[~np.isnan(np.unique((fire_prog_array)))])

        # test fire progression for the very first fire date
        fire_prog, cluster_with_correct_progression = self.fire_progression_v2(fire_prog_array)
        # save the pixel patches for which fire progression is normal 
        first_cluster = cluster_with_correct_progression
        
        # if there are more than two burn dates and the number of iterations is set to at least 2, 
        # a loop is needed to ckeck if the fire progresses normally over multiple dates
        if (data_differences_init_fire > 2) & (number_of_max_iterations>1):
            # check if the remaining number of fire dates is smaller than the number of iterations
            # since a lot of fires only last for two fire dates and the number of iterations is set to 3,
            # this happens often and it is therefore checked first 
            if ((data_differences_init_fire-2) < (number_of_max_iterations-1)):
                # check for fire progression for all the different dates 
                fire_prog, ori_cluster_list = self.fire_prog_iterations(first_cluster, fire_prog, cluster_with_correct_progression, 
                                                                (data_differences_init_fire-2))      
                return fire_prog, ori_cluster_list

            # if the remaining number of fire dates is larger than the number of iterations, start here
            else:
                # check for fire progression for specified number of iterations
                fire_prog, ori_cluster_list = self.fire_prog_iterations(first_cluster, fire_prog, cluster_with_correct_progression, 
                                                                (number_of_max_iterations-1))    
                return fire_prog, ori_cluster_list
        
        # if there are only two burn dates, return results of the initial fire progression check
        else:
            return fire_prog, first_cluster     