import xarray as xr
from datetime import datetime,timedelta
import os
import calendar

class iHARPExecuter():
    def __init__(self):
        self.current_directory = '/data/iHARPVFullStack/iHARPV/api/'
        self.dataDirectory = '/data/era5/'
        self.variable, self.variable_unit,self.variable_long_name,self.time_resolution,self.time_agg_method,self.spatial_resolution,self.spatial_agg_method = "","","","","","",""
        self.frames_directory = os.path.join(self.current_directory, 'assets/frames')
        self.min_lat,self.max_lat,self.min_lon,self.max_lon = -90,90,-180,180
    #This function is used in plotting
    def extractDateTimeInfo(self,start_date_time_str, end_date_time_str):
        # Parse start date-time string
        start_date_time = datetime.fromisoformat(start_date_time_str.rstrip("Z"))
        # Parse end date-time string
        end_date_time = datetime.fromisoformat(end_date_time_str.rstrip("Z"))

        # Extract information
        self.selected_month_name_start = str(start_date_time.strftime("%B")).lower()
        self.selected_month_name_end = str(end_date_time.strftime("%B")).lower()
        self.selected_month_start = datetime.strptime(self.selected_month_name_start, '%B').month
        self.selected_month_end = datetime.strptime(self.selected_month_name_end, '%B').month
        self.selected_hour_start = start_date_time.strftime("%H")
        self.selected_hour_end = end_date_time.strftime("%H")
        self.selected_day_start = start_date_time.strftime("%d").lstrip("0")
        self.selected_day_end = end_date_time.strftime("%d").lstrip("0")
        self.selected_year_start = start_date_time.strftime("%Y")
        self.selected_year_end = end_date_time.strftime("%Y")
        return    

    def getRaster(
        self,
        variable: str,
        start_datetime: str,
        end_datetime: str,
        time_resolution: str,  # e.g., "hour", "day", "month", "year"
        time_agg_method: str,  # e.g., "mean", "max", "min"
        south: float,
        north: float,
        west: float,
        east: float,
        spatial_resolution: float,  # e.g., 0.25, 0.5, 1.0, 2.5, 5.0
        spatial_agg_method: str,  # e.g., "mean", "max", "min"
    ) -> xr.Dataset:
        
        self.min_lat, self.max_lat, self.min_lon, self.max_lon,self.temporalLevel,self.time_agg_method,self.variable = float(south),float(north),float(west),float(east), time_resolution,time_agg_method,variable
        self.spatial_resolution,self.spatial_agg_method = float(spatial_resolution),spatial_agg_method
        if self.min_lon > self.max_lon:
            self.max_lon, self.min_lon = self.min_lon, self.max_lon
            self.max_lat, self.min_lat = self.min_lat, self.max_lat
        self.lon_range,self.lat_range = slice(self.min_lon, self.max_lon),slice(self.max_lat, self.min_lat) 
        self.startDateTime,self.endDateTime= start_datetime[:-11],end_datetime[:-11]
        start_date_counter = datetime.strptime(self.startDateTime, "%Y-%m-%dT%H")
        end_date_counter = datetime.strptime(self.endDateTime, "%Y-%m-%dT%H")
        self.selected_year_start,self.selected_year_end  = start_date_counter.strftime("%Y"),end_date_counter.strftime("%Y")
        print(self.selected_year_start,self.selected_year_end )
        if self.time_resolution=="Hourly":
            #CASE #1: Requested hours in more than more than one year #This is EXPENSIVE CASE
            if (str(self.selected_year_start) != str(self.selected_year_end)) :
                ds_list = []
                current_date = start_date_counter
                while current_date <= end_date_counter:
                    data_location = (self.dataDirectory+ 'raw/'+ self.variable+'/'+ self.variable+'-'+str(current_date.year)+'.nc')
                    # For the first year, retrieve data from selected start datetime till the end of the year
                    if current_date.year == self.selected_year_start:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(time=slice(self.startDateTime,None),longitude=self.lon_range, latitude=self.lat_range)
                    # For the last year, retrieve data from selected start datetime till the end of the year
                    elif current_date.year == self.selected_year_end:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(time=slice(None,self.endDateTime),longitude=self.lon_range, latitude=self.lat_range)
                    # For all years between get all hours in thos years
                    else:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(longitude=self.lon_range, latitude=self.lat_range)
                        ds_list.append(ds)
                    current_date += timedelta(years=1)
                ds = xr.concat(ds_list, dim="time")

            #CASE #2 Requested range of hours in one year
            else:
                data_location = (self.dataDirectory+'raw/'+  self.variable+'/'+ self.variable+'-'+str(self.selected_year_start)+'.nc')
                ds = xr.open_dataset(data_location, engine="netcdf4").sel(time=slice(self.startDateTime,self.endDateTime),longitude=self.lon_range, latitude=self.lat_range)
                
        else:
            print(self.dataDirectory+'preprocessed/' +self.variable+'/'+ 'combined_'+str(self.time_resolution).lower() + str(self.time_agg_method)+"_2014_2023.nc")
            data_location = self.dataDirectory+'preprocessed/' +self.variable+'/'+ 'combined_'+str(self.time_resolution).lower() + str(self.time_agg_method)+"_2014_2023.nc"
            ds = xr.open_dataset(data_location).sel(time=slice(self.startDateTime,self.endDateTime),longitude=self.lon_range, latitude=self.lat_range)
        
        if self.spatial_resolution!=0.25:
            coarsened = ds.coarsen(latitude=self.spatial_resolution/0.25,longitude=self.spatial_resolution/0.25)
            if self.spatial_agg_method=='mean':
                ds = coarsened.mean()
            elif self.spatial_agg_method=='min':
                ds = coarsened.min()
            elif self.spatial_agg_method=='max':
                ds = coarsened.max()    
        self.raster = ds

        return self.raster
        


    def viewTimeSeries(self, raster: xr.Dataset):

        pass


    def viewHeatMap(self, raster: xr.Dataset):
        pass


    def findArea(
        self,
        raster: xr.Dataset,
        agg_method: str,  # e.g., "mean", "max", "min"
        predicate: str,  # e.g., ">", "<", "==", "!=", ">=", "<="
        value: float,
    ) -> xr.Dataset:
        """
        add a new variable spatial_mask(lat, lon) to the dataset
        E.g., 
        dimensions: time, lat, lon
        variable: t2m(time, lat, lon)      (24, 721, 1440) --agg-> (721, 1440)
                  spatial_mask(lat, lon)   (721, 1440): boolean mask
        """
        pass


    def findTime(
        self,
        raster: xr.Dataset,
        agg_method: str,  # e.g., "mean", "max", "min"
        predicate: str,  # e.g., ">", "<", "==", "!=", ">=", "<="
        value: float,
    ) -> xr.Dataset:
        """
        add a new variable time_mask(time) to the dataset
        E.g., 
        dimensions: time, lat, lon
        variable: t2m(time, lat, lon)   (24, 721, 1440) --agg-> (24)
                  temporal_mask(time)   (24) : (true, true, false, ...)

        """
        pass

