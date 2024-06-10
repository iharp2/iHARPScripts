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
        time_resolution: str,  # e.g., "hourly", "daily", "monthly", "yearly"
        time_agg_method: str,  # e.g., "mean", "max", "min"
        south: float,
        north: float,
        west: float,
        east: float,
        spatial_resolution: float,  # e.g., 0.25, 0.5, 1.0, 2.5, 5.0
    ):
        if variable == "2m Temperature":
            variable = "t2m"
        elif variable == "Surface Pressure":
            variable = "sp"
        elif variable == "Total Precipitation":
            variable = "tp"
        self.variable = str(variable).lower()
        self.startDateTime, self.endDateTime = start_datetime[:-11], end_datetime[:-11]
        self.extract_date_time_info(self.startDateTime, self.endDateTime)
        self.time_resolution, self.time_agg_method = (
            str(time_resolution).lower(),
            str(time_agg_method).lower(),
        )
        self.min_lat = float(south)
        self.max_lat = float(north)
        self.min_lon = float(west)
        self.max_lon = float(east)
        self.spatial_resolution = float(spatial_resolution)
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
        self.variable_unit = raster[self.variable].attrs.get(
            "units", "No unit attribute found"
        )
        self.variable_long_name = raster[self.variable].attrs.get(
            "long_name", "No long name attribute found"
        )
        temperature_data = None
        image_location = os.path.join(
            self.current_directory, "assets/timeSeriesResult.json"
        )
        time_units, temps = [], []
        for time_step in raster.time:
            date = time_step.values
            dateTime = str(date)[:13]
            # Extract month and day from each time step
            year = time_step.values.astype("M8[Y]").astype("O").year
            month = time_step.values.astype("M8[M]").astype("O").month
            month_name = calendar.month_name[month]
            day = time_step.values.astype("M8[D]").astype("O").day
            if self.time_resolution == "hourly":
                time_unit = dateTime
            elif self.time_resolution == "daily":
                time_unit = str(year) + "-" + str(month) + "-" + str(day)
            elif self.time_resolution == "monthly":
                time_unit = str(year) + "-" + str(month)
            elif self.time_resolution == "yearly":
                time_unit = str(year)
            if time_unit in time_units:
                time_units.append(time_unit + "*")
            else:
                time_units.append(time_unit)

            xArray = raster[self.variable].sel(time=time_step)
            temperature_data = np.array(xArray.values)
            if self.time_agg_method == "min":
                # print("Entered Minimum raw")
                temp = np.min(temperature_data)
            elif self.time_agg_method == "max":
                # print("Entered Maximum raw")
                temp = np.max(temperature_data)
            else:
                temp = np.average(temperature_data)

            temps.append(temp)

        ytitle = str(self.variable_long_name + " [" + self.variable_unit + "]")
        # Naming the xtitle and title of the figure depends on teh time_resolution
        if self.time_resolution == "hourly":
            xtitle = self.time_resolution
            title = (
                self.time_agg_method
                + " "
                + str(self.variable_long_name)
                + " Variation "
                + self.time_resolution
                + " in "
                + str(month_name)
                + "-"
                + str(day)
                + "-"
                + str(year)
            )

        elif self.time_resolution == "daily":
            xtitle = "Year-Month-Day"
            if self.selected_month_name_start != self.selected_month_name_end:
                month_name = (
                    "("
                    + self.selected_month_name_start
                    + "~"
                    + self.selected_month_name_end
                    + ")"
                )
            title = (
                str(self.variable_long_name)
                + " Variation Daily in "
                + str(month_name)
                + "-"
                + str(year)
            )
        elif self.time_resolution == "monthly":
            xtitle = "Year-Month"
            if self.selected_year_start != self.selected_year_end:
                year = (
                    "("
                    + str(self.selected_year_start)
                    + "~"
                    + str(self.selected_year_end)
                    + ")"
                )
            else:
                year = self.selected_year_start
            title = str(self.variable_long_name) + " Variation Monthly in " + str(year)
        elif self.time_resolution == "yearly":
            xtitle = "Year"
            if self.selected_year_start != self.selected_year_end:
                year = (
                    "("
                    + str(self.selected_year_start)
                    + "~"
                    + str(self.selected_year_end)
                    + ")"
                )
            else:
                year = self.selected_year_start
            title = str(self.variable_long_name) + " Variation Yearly in " + str(year)

        trace = go.Scatter(
            x=time_units,
            y=temps,
            mode="lines+markers",
            name=self.time_agg_method + " " + self.variable_long_name,
            marker=dict(symbol="circle"),
        )
        fig = go.Figure(data=[trace])
        fig.update_xaxes(
            mirror=True,
            ticks="outside",
            showline=True,
            linecolor="black",
            gridcolor="lightgrey",
        )
        fig.update_yaxes(
            mirror=True,
            ticks="outside",
            showline=True,
            linecolor="black",
            gridcolor="lightgrey",
        )
        fig.update_layout(
            plot_bgcolor="white",
            xaxis_tickformat="%m-%d H:%H",
            # xaxis_rangeslider_visible=True,
            title=dict(
                text=title,
                font=dict(size=15, color="black", weight="bold"),
                x=0.5,
            ),
            xaxis_title=dict(text=xtitle, font=dict(size=12)),
            yaxis_title=dict(text=ytitle, font=dict(size=12)),
            modebar_remove=[
                # "toImage",
                "pan",
                "zoomIn",
                "zoomOut",
                "resetViewMapbox",
                "lasso2d",
                "select",
                "zoom",
                "autoScale",
                # "resetScale",
            ],
            margin={"r": 10, "t": 40, "b": 15, "l": 50},
            showlegend=True,
            legend=dict(
                font=dict(size=7),
                xanchor="right",  # Anchors the legend's x position
                yanchor="bottom",  # Anchors the legend's y position
            ),
            height=300,  # Adjust the height of the graph
            width=750,
        )

        # Update x-axis ticks
        # fig.update_xaxes(tickangle=10)
        fig_json = fig.to_json()
        with open(image_location, "w") as f:
            f.write(fig_json)
        # Read the JSON file
        with open(image_location, "r") as json_file:
            data = json.load(json_file)
        return data

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
