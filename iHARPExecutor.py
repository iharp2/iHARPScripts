import os
from datetime import datetime, timedelta
import calendar
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import plotly.graph_objs as go
import plotly.io as pio
import json
import numpy as np
import xarray as xr
import json
import pandas as pd
import matplotlib.pyplot as plt
from django.http import JsonResponse
from PIL import Image
import base64
from io import BytesIO
from matplotlib.animation import FuncAnimation
import plotly.express as px
import geopandas as gpd
from shapely.geometry import Polygon
from rest_framework.response import Response
from datetime import datetime, timedelta
import os
import numpy as np
import calendar
import matplotlib.dates as mdates
import plotly.graph_objs as go
import plotly.io as pio
from dateutil.relativedelta import relativedelta
from base64 import b64encode

try:
    from StringIO import StringIO

    py3 = False
except ImportError:
    from io import StringIO, BytesIO

    py3 = True
from random import random


class iHARPExecuter:
    def __init__(self):
        self.current_directory = "/data/iHARPVFullStack/iHARPV/api/"
        self.dataDirectory = "/data/era5/"
        (
            self.variable,
            self.variable_unit,
            self.variable_long_name,
            self.time_resolution,
            self.time_agg_method,
            self.spatial_resolution,
            self.spatial_agg_method,
        ) = ("", "", "", "", "", "", "")
        self.frames_directory = os.path.join(self.current_directory, "assets/frames")
        self.min_lat, self.max_lat, self.min_lon, self.max_lon = -90, 90, -180, 180
        var = random()
        print(var)

    # This function is used in plotting
    def extract_date_time_info(self, start_date_time_str, end_date_time_str):
        # Parse start date-time string
        start_date_time = datetime.fromisoformat(start_date_time_str.rstrip("Z"))
        # Parse end date-time string
        end_date_time = datetime.fromisoformat(end_date_time_str.rstrip("Z"))

        # Extract information
        self.selected_month_name_start = str(start_date_time.strftime("%B")).lower()
        self.selected_month_name_end = str(end_date_time.strftime("%B")).lower()
        self.selected_month_start = datetime.strptime(
            self.selected_month_name_start, "%B"
        ).month
        self.selected_month_end = datetime.strptime(
            self.selected_month_name_end, "%B"
        ).month
        self.selected_hour_start = start_date_time.strftime("%H")
        self.selected_hour_end = end_date_time.strftime("%H")
        self.selected_day_start = start_date_time.strftime("%d").lstrip("0")
        self.selected_day_end = end_date_time.strftime("%d").lstrip("0")
        self.selected_year_start = start_date_time.strftime("%Y")
        self.selected_year_end = end_date_time.strftime("%Y")
        return

    def getRaster_old(
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
    ):
        if variable == "2m Temperature":
            variable = "t2m"
        elif variable == "Surface Pressure":
            variable = "sp"
        elif variable == "Total Precipitation":
            variable = "tp"
        (
            self.min_lat,
            self.max_lat,
            self.min_lon,
            self.max_lon,
            self.time_resolution,
            self.time_agg_method,
            self.variable,
        ) = (
            float(south),
            float(north),
            float(west),
            float(east),
            time_resolution,
            time_agg_method,
            variable,
        )
        self.spatial_resolution, self.spatial_agg_method = (
            float(spatial_resolution),
            spatial_agg_method,
        )
        if self.min_lon > self.max_lon:
            self.max_lon, self.min_lon = self.min_lon, self.max_lon
            self.max_lat, self.min_lat = self.min_lat, self.max_lat
        self.lon_range, self.lat_range = slice(self.min_lon, self.max_lon), slice(
            self.max_lat, self.min_lat
        )
        self.startDateTime, self.endDateTime = start_datetime[:-11], end_datetime[:-11]
        self.extract_date_time_info(self.startDateTime, self.endDateTime)
        start_date_counter = datetime.strptime(self.startDateTime, "%Y-%m-%dT%H")
        end_date_counter = datetime.strptime(self.endDateTime, "%Y-%m-%dT%H")
        if self.time_resolution == "Hourly":
            # CASE #1: Requested hours in more than more than one year #This is EXPENSIVE CASE
            if str(self.selected_year_start) != str(self.selected_year_end):
                ds_list = []
                current_date = start_date_counter
                while current_date <= end_date_counter:
                    data_location = (
                        self.dataDirectory
                        + "raw/"
                        + self.variable
                        + "/"
                        + self.variable
                        + "-"
                        + str(current_date.year)
                        + ".nc"
                    )
                    # For the first year, retrieve data from selected start datetime till the end of the year
                    if current_date.year == self.selected_year_start:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            time=slice(self.startDateTime, None),
                            longitude=self.lon_range,
                            latitude=self.lat_range,
                        )
                    # For the last year, retrieve data from selected start datetime till the end of the year
                    elif current_date.year == self.selected_year_end:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            time=slice(None, self.endDateTime),
                            longitude=self.lon_range,
                            latitude=self.lat_range,
                        )
                    # For all years between get all hours in thos years
                    else:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            longitude=self.lon_range, latitude=self.lat_range
                        )
                        ds_list.append(ds)
                    current_date += timedelta(years=1)
                ds = xr.concat(ds_list, dim="time")

            # CASE #2 Requested range of hours in one year
            else:
                data_location = (
                    self.dataDirectory
                    + "raw/"
                    + self.variable
                    + "/"
                    + self.variable
                    + "-"
                    + str(self.selected_year_start)
                    + ".nc"
                )
                ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )

        else:

            data_location = (
                self.dataDirectory
                + "preprocessed/"
                + self.variable
                + "/"
                + "combined_"
                + str(self.time_resolution).lower()
                + "_"
                + str(self.time_agg_method)
                + "_2014_2023.nc"
            )
            ds = xr.open_dataset(data_location).sel(
                time=slice(self.startDateTime, self.endDateTime),
                longitude=self.lon_range,
                latitude=self.lat_range,
            )

        if self.spatial_resolution != 0.25:
            coarsened = ds.coarsen(
                latitude=int(self.spatial_resolution / 0.25),
                longitude=int(self.spatial_resolution / 0.25),
                boundary="trim",
            )
            if self.spatial_agg_method == "mean":
                ds = coarsened.mean()
            elif self.spatial_agg_method == "min":
                ds = coarsened.min()
            elif self.spatial_agg_method == "max":
                ds = coarsened.max()
        self.raster = ds

        return self.raster

    def checkDataSupportance(self):
        if self.variable == "t2m":
            return True
        if self.variable == "tp":
            if (
                self.time_resolution == "hourly"
                and int(self.selected_year_start) < 2019
            ):
                return False
            return True
        if self.variable == "sp":
            if self.time_resolution == "hourly":
                return False
            else:  # time_resolution == "daily", "monthly", "yearly"
                if self.min_lat >= 30:  # part 1: lat [90, 30]
                    return True
                elif self.max_lat <= 30:  # part 2: lat [30, -90]
                    if float(self.spatial_resolution) >= 1:
                        return True
                    else:
                        return False
                else:  # Deny inter-spatial-partition query for now
                    return False

    def genFileList(self):
        if self.variable == "t2m" or self.variable == "tp":
            if self.time_resolution == "hourly":
                return [
                    f"/data/era5/raw/{self.variable}/{self.variable}-{year}.nc"
                    for year in range(
                        int(self.selected_year_start), int(self.selected_year_end) + 1
                    )
                ]
            else:
                return [
                    f"/data/era5/preprocessed/{self.variable}/combined_{self.time_resolution}_{self.time_agg_method}_2014_2023.nc"
                ]
        if self.variable == "sp":
            if self.min_lat >= 30:  # part 1: lat [90, 30]
                return [
                    f"/data/era5/preprocessed/{self.variable}/combined_{self.time_resolution}_{self.time_agg_method}_2014_2023_part1.nc"
                ]
            elif self.max_lat <= 30:  # part 2: lat [30, -90]
                return [
                    f"/data/era5/preprocessed/{self.variable}/combined_{self.time_resolution}_{self.time_agg_method}_2014_2023_part2.nc"
                ]
            else:
                raise ValueError(
                    "iHARPV: Inter-spatial-partition query is not supported yet."
                )

    def spatialDownsample(self, ds):
        return ds.sel(
            latitude=np.arange(self.max_lat, self.min_lat, -self.spatial_resolution),
            longitude=np.arange(self.min_lon, self.max_lon, self.spatial_resolution),
            method="nearest",
        )

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
    ):
        self.variable = variable
        self.startDateTime, self.endDateTime = start_datetime[:-11], end_datetime[:-11]
        self.extract_date_time_info(self.startDateTime, self.endDateTime)
        self.time_resolution, self.time_agg_method = time_resolution, time_agg_method
        self.min_lat = south
        self.max_lat = north
        self.min_lon = west
        self.max_lon = east
        self.spatial_resolution = spatial_resolution

        if not self.checkDataSupportance():
            raise ValueError("iHARPV: Query range or resolution not supported")

        file_list = self.genFileList()
        ds_list = []
        for file in file_list:
            ds = xr.open_dataset(file, engine="netcdf4").sel(
                time=slice(self.startDateTime, self.endDateTime),
                latitude=slice(self.max_lat, self.min_lat),
                longitude=slice(self.min_lon, self.max_lon),
            )
            ds_list.append(ds)
        ds = xr.concat(ds_list, dim="time")

        if self.spatial_resolution != 0.25:
            ds = self.spatialDownsample(ds)

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

    def getHeatMapScript(
        self,
        variable,
        startDateTime,
        endDateTime,
        temporalLevel,
        aggLevel,
        north,
        east,
        south,
        west,
    ):
        self.min_lat, self.max_lat, self.min_lon, self.max_lon = (
            float(south),
            float(north),
            float(west),
            float(east),
        )
        self.temporalLevel = temporalLevel
        self.variable = variable
        self.aggLevel = aggLevel
        if self.min_lon > self.max_lon:
            self.max_lon, self.min_lon = self.min_lon, self.max_lon
            self.max_lat, self.min_lat = self.min_lat, self.max_lat

        # Specify the longitude and latitudenrange
        self.lon_range = slice(self.min_lon, self.max_lon)
        self.lat_range = slice(self.max_lat, self.min_lat)
        self.startDateTime = startDateTime[:-11]
        self.endDateTime = endDateTime[:-11]
        self.extract_date_time_info(self.startDateTime, self.endDateTime)
        ## Loading Initial Data
        videoPath = os.path.join(self.current_directory, "assets/heatmapVideo.mp4")
        # print(videoPath)
        if variable == "2m Temperature":
            self.variable = "t2m"
        elif variable == "Surface Pressure":
            self.variable = "sp"
        elif variable == "Total Precipitation":
            self.variable = "tp"
        if temporalLevel == "Hourly":
            # CASE #1: Requested more than more than one year of hours
            if str(self.selected_year_start) != str(self.selected_year_end):
                ds_list = []
                start_date = datetime.strptime(self.startDateTime, "%Y-%m-%dT%H")
                end_date = datetime.strptime(self.endDateTime, "%Y-%m-%dT%H")
                current_date = start_date
                while current_date <= end_date:
                    data_location = (
                        self.dataDirectory
                        + self.variable
                        + "/"
                        + self.variable
                        + "-"
                        + str(current_date.year)
                        + ".nc"
                    )
                    # print(data_location)
                    # For the first year, retrieve data from selected start datetime till the end of the year
                    if current_date.year == self.selected_year_start:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            time=slice(self.startDateTime, None),
                            longitude=self.lon_range,
                            latitude=self.lat_range,
                        )
                    # For the last year, retrieve data from selected start datetime till the end of the year
                    elif current_date.year == self.selected_year_end:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            time=slice(None, self.endDateTime),
                            longitude=self.lon_range,
                            latitude=self.lat_range,
                        )
                    # For all years between get all hours in thos years
                    else:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            longitude=self.lon_range, latitude=self.lat_range
                        )
                    ds_list.append(ds)
                    current_date += timedelta(years=1)
                self.ds_daily = xr.concat(ds_list, dim="time")

            # CASE #2 Requested range of hours in one day
            else:
                data_location = (
                    self.dataDirectory
                    + "raw/"
                    + self.variable
                    + "/"
                    + self.variable
                    + "-"
                    + str(self.selected_year_start)
                    + ".nc"
                )
                self.ds_daily = xr.open_dataset(data_location, engine="netcdf4").sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )

            # update_progress(100,'Retrieving Time Series..Finished:')
        else:
            start_date = datetime.strptime(self.startDateTime, "%Y-%m-%dT%H")
            end_date = datetime.strptime(self.endDateTime, "%Y-%m-%dT%H")
            if temporalLevel == "Daily":
                data_location = (
                    self.dataDirectory
                    + "preprocessed/"
                    + self.variable
                    + "/"
                    + "combined_daily_"
                )
                # self.endDateTime = end_date.replace(day=end_date.day + 1)
            elif temporalLevel == "Monthly":
                data_location = (
                    self.dataDirectory
                    + "preprocessed/"
                    + self.variable
                    + "/"
                    + "combined_monthly_"
                )
                # self.endDateTime = end_date.replace(month=end_date.month + 1)
            elif temporalLevel == "Yearly":
                data_location = (
                    self.dataDirectory
                    + "preprocessed/"
                    + self.variable
                    + "/"
                    + "combined_yearly_"
                )
                # self.endDateTime = end_date.replace(year=end_date.year + 1)
            # Getting Maximum Data Per Month
            # self.endDateTime
            if self.aggLevel == "Maximum":
                self.daily_temp_max = xr.open_dataset(
                    data_location + "max_2014_2023.nc"
                ).sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )

            elif self.aggLevel == "Minimum":
                self.daily_temp_min = xr.open_dataset(
                    data_location + "min_2014_2023.nc"
                ).sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )

            elif self.aggLevel == "Mean":
                self.daily_temp_avg = xr.open_dataset(
                    data_location + "mean_2014_2023.nc"
                ).sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )
            # Getting Minimimum Data Per Month
            # Getting Average Data Per Month

        self.HeatMapTemps = []
        # print("Building HeatMap Hourly")
        # print(self.frames_directory)
        for file_name in os.listdir(self.frames_directory):
            # print(file_name)
            if file_name.endswith(".png"):
                file_path_full = os.path.join(self.frames_directory, file_name)
                # print(file_path_full)
                os.remove(file_path_full)
        # print(self.ds_daily.time)
        if self.temporalLevel == "Hourly":
            for time_step in self.ds_daily.time:
                xArray = self.ds_daily[self.variable].sel(time=time_step)
                values = xArray.values
                xArrayCopy = xArray
                xArrayCopy.values = values
                self.HeatMapTemps.append(xArrayCopy)
        else:
            if self.aggLevel == "Mean":
                daily_temp_heatmap = self.daily_temp_avg[self.variable]
            elif self.aggLevel == "Minimum":
                daily_temp_heatmap = self.daily_temp_min[self.variable]
            elif self.aggLevel == "Maximum":
                daily_temp_heatmap = self.daily_temp_max[self.variable]
            for time_Step in daily_temp_heatmap.time:
                daily = daily_temp_heatmap.sel(time=time_Step)
                values = daily.values
                xArrayCopy = daily
                xArrayCopy.values = values
                self.HeatMapTemps.append(xArrayCopy)
        # TODO: We need to figure out how to update the progress bar
        # update_progress(75,'Retrieving Time Series..Progress:')
        counter = 0
        aggregationMethod = self.aggLevel
        total = len(self.HeatMapTemps)
        index = -1
        self.variable_unit = self.HeatMapTemps[0].attrs.get(
            "units", "No unit attribute found"
        )
        self.variable_long_name = self.HeatMapTemps[0].attrs.get(
            "long_name", "No long name attribute found"
        )
        for data in self.HeatMapTemps:
            # print(data)
            progress_percentage = (index + 1) / total * 100
            # print(f"Progress: {progress_percentage:.2f}%")
            # update_progress(int(progress_percentage),'Building HeatMap..Progress:')
            index += 1
            time_step = data.time
            # data.attrs["GRIB_units"] = "C"
            # data.attrs["units"] = "C"
            hour = time_step.dt.hour.item()
            hour_am_pm = datetime.strptime(str(hour), "%H").strftime("%I %p")

            month = time_step.values.astype("M8[M]").astype("O").month
            month_name = calendar.month_name[month]
            day = time_step.values.astype("M8[D]").astype("O").day
            year = time_step.dt.year.item()
            fig, ax = plt.subplots(figsize=(12, 9))
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            if self.temporalLevel == "Hourly":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(month_name)}  {str(day)} at : {hour_am_pm}",
                    ax=ax,
                )
                ax.set_title(
                    f"Hourly {self.variable_long_name} Variation - {str(month_name)+ str(day)} at : {hour_am_pm}"
                )

            elif self.temporalLevel == "Daily":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(month_name)+str(day)} at : {hour_am_pm}",
                    ax=ax,
                )
                ax.set_title(
                    f"Daily {aggregationMethod} {self.variable_long_name} Variation - on {str(month_name)} - {str(day)}"
                )
            elif self.temporalLevel == "Monthly":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(month_name)+str(year)}",
                    ax=ax,
                )
                ax.set_title(
                    f"Monthly {aggregationMethod} {self.variable_long_name} Variation - on {str(month_name)} - {str(year)}"
                )
            elif self.temporalLevel == "Yearly":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(year)}", ax=ax
                )
                ax.set_title(
                    f"Yearly {aggregationMethod} {self.variable_long_name} Variation - on {str(year)}"
                )

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.tick_params(
                axis="both",
                which="both",
                bottom=False,
                top=False,
                left=False,
                right=False,
            )
            ax.grid(True, linestyle="--", alpha=0.2)  # Add grid lines

            # Adjust margins
            plt.subplots_adjust(
                left=0.1, right=0.9, top=0.9, bottom=0.1
            )  # Adjust as needed

            # Save the plo t
            output_file = self.frames_directory + f"/{counter}.png"
            counter += 1
            plt.savefig(output_file)
            plt.close(fig)  # Close the figure to free up memory
        # Create a video file from the images
        os.system(
            f'ffmpeg -framerate 2 -start_number 0 -i  "{self.frames_directory}/%d.png" -c:v libx264 -r 30 "{videoPath}" -y'
        )

        # Open the image file
        with open(videoPath, "rb") as f:
            video_data = f.read()

        video_file_bytes = BytesIO(video_data)
        # After you send the video go and delete all frames
        for file_name in os.listdir(self.frames_directory):
            # print(file_name)
            if file_name.endswith(".png"):
                file_path_full = os.path.join(self.frames_directory, file_name)
                # print(file_path_full)
                os.remove(file_path_full)
        # Return the response as JSON
        return video_file_bytes
        # return 1

    def getAreas(
        self,
        variable,
        startDateTime,
        endDateTime,
        temporalLevel,
        aggLevel,
        north,
        east,
        south,
        west,
    ):
        self.min_lat, self.max_lat, self.min_lon, self.max_lon = (
            float(south),
            float(north),
            float(west),
            float(east),
        )
        self.temporalLevel = temporalLevel
        self.variable = variable
        self.aggLevel = aggLevel
        if self.min_lon > self.max_lon:
            self.max_lon, self.min_lon = self.min_lon, self.max_lon
            self.max_lat, self.min_lat = self.min_lat, self.max_lat

        # Specify the longitude and latitudenrange
        self.lon_range = slice(self.min_lon, self.max_lon)
        self.lat_range = slice(self.max_lat, self.min_lat)
        self.startDateTime = startDateTime[:-11]
        self.endDateTime = endDateTime[:-11]
        self.extract_date_time_info(self.startDateTime, self.endDateTime)
        json_file_path = os.path.join(self.current_directory, "assets/areasData.json")
        if variable == "2m Temperature":
            self.variable = "t2m"
        elif variable == "Surface Pressure":
            self.variable = "sp"
        elif variable == "Total Precipitation":
            self.variable = "tp"
        if temporalLevel == "Hourly":
            # CASE #1: Requested more than more than one year of hours
            if str(self.selected_year_start) != str(self.selected_year_end):
                ds_list = []
                start_date = datetime.strptime(self.startDateTime, "%Y-%m-%dT%H")
                end_date = datetime.strptime(self.endDateTime, "%Y-%m-%dT%H")
                current_date = start_date
                while current_date <= end_date:
                    data_location = (
                        self.dataDirectory
                        + self.variable
                        + "/"
                        + self.variable
                        + "-"
                        + str(current_date.year)
                        + ".nc"
                    )
                    # print(data_location)
                    # For the first year, retrieve data from selected start datetime till the end of the year
                    if current_date.year == self.selected_year_start:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            time=slice(self.startDateTime, None),
                            longitude=self.lon_range,
                            latitude=self.lat_range,
                        )
                    # For the last year, retrieve data from selected start datetime till the end of the year
                    elif current_date.year == self.selected_year_end:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            time=slice(None, self.endDateTime),
                            longitude=self.lon_range,
                            latitude=self.lat_range,
                        )
                    # For all years between get all hours in thos years
                    else:
                        ds = xr.open_dataset(data_location, engine="netcdf4").sel(
                            longitude=self.lon_range, latitude=self.lat_range
                        )
                    ds_list.append(ds)
                    current_date += timedelta(years=1)
                self.ds_daily = xr.concat(ds_list, dim="time")

            # CASE #2 Requested range of hours in one day
            else:
                data_location = (
                    self.dataDirectory
                    + "raw/"
                    + self.variable
                    + "/"
                    + self.variable
                    + "-"
                    + str(self.selected_year_start)
                    + ".nc"
                )
                self.ds_daily = xr.open_dataset(data_location, engine="netcdf4").sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )

            # update_progress(100,'Retrieving Time Series..Finished:')
        else:
            start_date = datetime.strptime(self.startDateTime, "%Y-%m-%dT%H")
            end_date = datetime.strptime(self.endDateTime, "%Y-%m-%dT%H")
            if temporalLevel == "Daily":
                data_location = (
                    self.dataDirectory
                    + "preprocessed/"
                    + self.variable
                    + "/"
                    + "combined_daily_"
                )
                # self.endDateTime = end_date.replace(day=end_date.day + 1)
            elif temporalLevel == "Monthly":
                data_location = (
                    self.dataDirectory
                    + "preprocessed/"
                    + self.variable
                    + "/"
                    + "combined_monthly_"
                )
                # self.endDateTime = end_date.replace(month=end_date.month + 1)
            elif temporalLevel == "Yearly":
                data_location = (
                    self.dataDirectory
                    + "preprocessed/"
                    + self.variable
                    + "/"
                    + "combined_yearly_"
                )
                # self.endDateTime = end_date.replace(year=end_date.year + 1)
            if self.aggLevel == "Maximum":
                self.daily_temp_max = xr.open_dataset(
                    data_location + "max_2014_2023.nc"
                ).sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )
                self.ds_daily = self.daily_temp_max

            elif self.aggLevel == "Minimum":
                self.daily_temp_min = xr.open_dataset(
                    data_location + "min_2014_2023.nc"
                ).sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )
                self.ds_daily = self.daily_temp_min

            elif self.aggLevel == "Mean":
                self.daily_temp_avg = xr.open_dataset(
                    data_location + "mean_2014_2023.nc"
                ).sel(
                    time=slice(self.startDateTime, self.endDateTime),
                    longitude=self.lon_range,
                    latitude=self.lat_range,
                )
                self.ds_daily = self.daily_temp_avg

        # Assume the surface pressure variable is named 'sp'
        self.ds_daily = self.ds_daily.coarsen(
            latitude=8, longitude=8, boundary="trim"
        ).mean()

        surface_pressure = self.ds_daily[self.variable]
        print("Total number of values:", surface_pressure.size)

        # Define the threshold for surface pressure
        threshold_max = 92000
        threshold_min = 90000

        # Create a boolean mask indicating whether each value meets the condition
        condition_met = (surface_pressure > threshold_min) & (
            surface_pressure < threshold_max
        )

        # Convert the DataArray and mask to DataFrames
        surface_pressure_df = surface_pressure.to_dataframe().reset_index()
        condition_met_df = condition_met.to_dataframe(
            name="condition_met"
        ).reset_index()

        # Merge the two DataFrames
        filtered_df = pd.merge(
            surface_pressure_df, condition_met_df, on=["time", "latitude", "longitude"]
        )
        df_sliced_json = filtered_df.head(1000).to_dict("records")
        print(filtered_df.head())
        # Print the total number of values and number of non-NaN values
        print("Total number of values in DataFrame:", filtered_df["sp"].size)
        print("Number of non-NaN values in DataFrame:", filtered_df["sp"].count())

        # Ensure that 'latitude' and 'longitude' are the coordinate names
        filtered_df = filtered_df.rename(
            columns={"lat": "latitude", "lon": "longitude"}
        )

        # Print the DataFrame to verify the results
        df = filtered_df
        df["latitude"] = df["latitude"] - 1
        df["longitude"] = df["longitude"] - 1
        df["latitude2"] = df["latitude"] + 2
        df["longitude2"] = df["longitude"] + 2
        gdf = gpd.GeoDataFrame(
            df,
            geometry=[
                Polygon([(x, y), (x, y2), (x2, y2), (x2, y)])
                for x, y, x2, y2 in zip(
                    df["longitude"], df["latitude"], df["longitude2"], df["latitude2"]
                )
            ],
        )
        # print("I reached here")
        fig = px.choropleth_mapbox(
            gdf,
            color_continuous_scale="Viridis",
            geojson=gdf.geometry,
            locations=gdf.index,
            color="condition_met",
            mapbox_style="white-bg",
            center={"lat": gdf["latitude"].mean(), "lon": gdf["longitude"].mean()},
            opacity=0.5,
            zoom=1,
        )

        fig.update_layout(
            mapbox_style="white-bg",
            mapbox_layers=[
                {
                    "below": "traces",
                    "sourcetype": "raster",
                    "sourceattribution": "United States Geological Survey",
                    "source": [
                        "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
                    ],
                }
            ],
            modebar_remove=[
                "toImage",
                "pan",
                "zoomIn",
                "zoomOut",
                "resetViewMapbox",
                "lasso2d",
                "select",
            ],
            # displaylogo=False,  # Hide the Plotly logo
            margin={"r": 0, "t": 0, "l": 0, "b": 0},
            showlegend=True,
            legend=dict(
                font=dict(size=7),
                x=1,  # Adjust the x position (0 to 1, 1 is far right)
                y=1,  # Adjust the y position (0 to 1, 1 is top)
                xanchor="right",  # Anchors the legend's x position
                yanchor="top",  # Anchors the legend's y position
            ),
        )

        fig_json = fig.to_json()
        with open(json_file_path, "w") as f:
            f.write(fig_json)
        # Read the JSON file
        with open(json_file_path, "r") as json_file:
            data = json.load(json_file)
        # Combine the two JSON objects into one dictionary
        # df_sliced_json = filtered_df.head(100).to_json(orient="records")

        combined_data = {"plotlyData": data, "dfData": df_sliced_json}
        # print(df_sliced_json)
        return combined_data
