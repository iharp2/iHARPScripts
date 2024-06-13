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
import operator

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
        time_units, temps_min, temps_max, temps_mean = [], [], [], []
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
            # if self.time_agg_method == "min":
            # print("Entered Minimum raw")
            temp_min = np.min(temperature_data)
            temps_min.append(temp_min)
            # elif self.time_agg_method == "max":
            # print("Entered Maximum raw")
            temp_max = np.max(temperature_data)
            temps_max.append(temp_max)
            # else:
            temp_mean = np.average(temperature_data)
            temps_mean.append(temp_mean)

        ytitle = str(self.variable_long_name + " [" + self.variable_unit + "]")
        # Naming the xtitle and title of the figure depends on teh time_resolution
        if self.time_resolution == "hourly":
            if self.selected_day_start != self.selected_hour_end:
                day = (
                    "("
                    + str(self.selected_day_start)
                    + "~"
                    + str(self.selected_day_end)
                    + ")"
                )
            title = (
                str(self.variable_long_name)
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
            # xtitle = "Year"
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

        trace1 = go.Scatter(
            x=time_units,
            y=temps_min,
            mode="lines+markers",
            name="min",
            marker=dict(symbol="circle"),
        )
        trace2 = go.Scatter(
            x=time_units,
            y=temps_max,
            mode="lines+markers",
            name="max",
            marker=dict(symbol="square"),
        )
        trace3 = go.Scatter(
            x=time_units,
            y=temps_mean,
            mode="lines+markers",
            name="mean",
            marker=dict(symbol="triangle-up"),
        )
        fig = go.Figure(data=[trace1, trace2, trace3])
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
            xaxis_tickformat="%Y-%m-%d<br>%H:00",
            # xaxis_rangeslider_visible=True,
            title=dict(
                text=title,
                font=dict(size=15, color="black", weight="bold"),
                x=0.5,
            ),
            xaxis_title=dict(text="Time", font=dict(size=12)),
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
                font=dict(size=14),
                xanchor="right",  # Anchors the legend's x position
                yanchor="bottom",  # Anchors the legend's y position
                x=1.14,
                y=0.7,
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
        videoPath = os.path.join(self.current_directory, "assets/heatmapVideo.mp4")
        self.HeatMapTemps = []
        for file_name in os.listdir(self.frames_directory):
            # print(file_name)
            if file_name.endswith(".png"):
                file_path_full = os.path.join(self.frames_directory, file_name)
                # print(file_path_full)
                os.remove(file_path_full)
        for time_step in raster.time:
            self.HeatMapTemps.append(raster[self.variable].sel(time=time_step))
        counter = 0
        index = -1
        self.variable_unit = self.HeatMapTemps[0].attrs.get(
            "units", "No unit attribute found"
        )
        self.variable_long_name = self.HeatMapTemps[0].attrs.get(
            "long_name", "No long name attribute found"
        )
        for data in self.HeatMapTemps:
            index += 1
            time_step = data.time
            hour = time_step.dt.hour.item()
            hour_am_pm = datetime.strptime(str(hour), "%H").strftime("%I %p")

            month = time_step.values.astype("M8[M]").astype("O").month
            month_name = calendar.month_name[month]
            day = time_step.values.astype("M8[D]").astype("O").day
            year = time_step.dt.year.item()
            fig, ax = plt.subplots(figsize=(20, 14))
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            if self.time_resolution == "hourly":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(month_name)}  {str(day)} at : {hour_am_pm}",
                    ax=ax,
                )
                ax.set_title(
                    f"{self.time_resolution} {self.variable_long_name} Variation - {str(month_name)+ str(day)} at : {hour_am_pm}"
                )

            elif self.time_resolution == "daily":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(month_name)+str(day)} at : {hour_am_pm}",
                    ax=ax,
                )
                ax.set_title(
                    f"{self.time_resolution} {self.time_agg_method} {self.variable_long_name} Variation - on {str(month_name)} - {str(day)}"
                )
            elif self.time_resolution == "monthly":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(month_name)+str(year)}",
                    ax=ax,
                )
                ax.set_title(
                    f"{self.time_resolution} {self.time_agg_method} {self.variable_long_name} Variation - on {str(month_name)} - {str(year)}"
                )
            elif self.time_resolution == "yearly":
                data.plot.imshow(
                    label=f"{self.variable_long_name} of {str(year)}", ax=ax
                )
                ax.set_title(
                    f"{self.time_resolution} {self.time_agg_method} {self.variable_long_name} Variation - on {str(year)}"
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
                left=0.05, right=1.08, top=0.97, bottom=0.08
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

    def findArea(
        self,
        raster: xr.Dataset,
        agg_method: str,  # e.g., "mean", "max", "min"
        predicate: str,  # e.g., ">", "<", "==", "!=", ">=", "<="
        value: float,
    ) -> xr.Dataset:

        json_file_path = os.path.join(self.current_directory, "assets/areasData.json")
        # Assume the surface pressure variable is named 'sp'
        # Dictionary mapping operator strings to their corresponding functions
        operator_mapping = {
            ">": operator.gt,
            "<": operator.lt,
            "=": operator.eq,  # Note that in Python, '==' is used for equality check
            "<=": operator.le,
            ">=": operator.ge,
            "!=": operator.ne,
        }
        # Aggregate over time dimension
        op_func = operator_mapping[predicate]
        variable = raster[self.variable]
        if agg_method == "any":
            # Apply the condition across all raster tables
            condition_test = op_func(variable, value)
            condition_met = (
                raster.where(condition_test)[self.variable]
                .notnull()
                .any(dim="time")
                .astype(bool)
            )
            var2 = raster[self.variable].isel(time=0)

        else:
            condition_test = op_func(variable, value)
            condition_met = (
                raster.where(condition_test)[self.variable]
                .notnull()
                .all(dim="time")
                .astype(bool)
            )
            var2 = raster[self.variable].isel(time=0)

        variable_df = var2.to_dataframe().reset_index()
        condition_met_df = condition_met.to_dataframe(
            name="condition_met"
        ).reset_index()
        # Merge the two DataFrames
        filtered_df = pd.merge(
            variable_df, condition_met_df, on=["latitude", "longitude"]
        )
        # Just send to the user the first 1000 values if it's more than this
        # Ensure that 'latitude' and 'longitude' are the coordinate names
        filtered_df = filtered_df.rename(columns={self.variable: "variable"})
        df_sliced_json = filtered_df.head(2000).to_dict("records")
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
        color_mapping = {True: "blue", False: "red"}

        fig = px.choropleth_mapbox(
            gdf,
            color_continuous_scale="Viridis",
            geojson=gdf.geometry,
            locations=gdf.index,
            hover_data={"condition_met": True, "variable": True},
            color="condition_met",
            mapbox_style="white-bg",
            center={"lat": gdf["latitude"].mean(), "lon": gdf["longitude"].mean()},
            opacity=0.5,
            zoom=1,
            color_discrete_map=color_mapping,
        )
        # fig.update_traces(
        #     hovertemplate="<b>Condition Met:</b> %{customdata[0]}<br><b>Variable:</b> %{customdata[1]}",
        #     customdata=gdf[
        #         ["condition_met", "variable"]
        #     ].values,  # Specify columns for customdata
        # )

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
        fig = ""
        return combined_data

    def findTime(
        self,
        raster: xr.Dataset,
        agg_method: str,  # e.g., "mean", "max", "min"
        predicate: str,  # e.g., ">", "<", "==", "!=", ">=", "<="
        value: float,
    ) -> xr.Dataset:
        json_file_path = os.path.join(self.current_directory, "assets/timesData.json")
        # Assume the surface pressure variable is named 'sp'
        # Dictionary mapping operator strings to their corresponding functions
        operator_mapping = {
            ">": operator.gt,
            "<": operator.lt,
            "=": operator.eq,  # Note that in Python, '==' is used for equality check
            "<=": operator.le,
            ">=": operator.ge,
            "!=": operator.ne,
        }
        # Aggregate over time dimension
        op_func = operator_mapping[predicate]
        variable = raster[self.variable]
        mytime = raster.time.values
        if agg_method == "any":
            # Apply the condition across all raster tables
            condition_test = op_func(variable, value)
            condition_met = (
                raster.where(condition_test)[self.variable]
                .notnull()
                .any(dim=["latitude", "longitude"])
                .astype(bool)
            )
        else:
            condition_test = op_func(variable, value)
            condition_met = (
                raster.where(condition_test)[self.variable]
                .notnull()
                .all(dim=["latitude", "longitude"])
                .astype(bool)
            )
        column_names = ["time"]

        time_df = pd.DataFrame(mytime, columns=column_names)
        condition_met_df = condition_met.to_dataframe(
            name="condition_met"
        ).reset_index()
        filtered_df = pd.merge(time_df, condition_met_df, on=["time"])
        df_sliced_json = filtered_df.head(2000).to_dict("records")
        condition_met = filtered_df["condition_met"]
        fig = go.Figure()
        filtered_df["value_numeric"] = filtered_df["condition_met"].astype(int)
        filtered_df["color"] = filtered_df["condition_met"].apply(
            lambda x: "blue" if x else "red"
        )

        fig.add_trace(
            go.Scatter(
                x=filtered_df["time"],
                y=filtered_df["value_numeric"],
                mode="lines+markers",
                marker=dict(
                    color=filtered_df["color"],
                    size=10,  # Adjust marker size as needed
                    line=dict(width=0.5),  # Adjust marker border width as needed
                ),
                line=dict(color="rgba(0, 0, 0, 0.6)"),  # Line color
                name="Condition Met",
            )
        )

        fig.update_yaxes(
            tickvals=[0, 1], ticktext=["False", "True"], title="Condition", side="left"
        )

        # Add lines to connect points
        fig.update_traces(mode="lines+markers")

        # Update layout to include secondary y-axis
        fig.update_layout(
            plot_bgcolor="white",
            xaxis_tickformat="%Y-%m-%d<br>%H:00",
            # xaxis_rangeslider_visible=True,
            title=dict(
                text="Find Times Result",
                font=dict(size=15, color="black", weight="bold"),
                x=0.5,
            ),
            xaxis_title=dict(text="Time", font=dict(size=12)),
            legend=dict(
                x=0,  # Position legend at the top-left corner
                y=1.35,
                bordercolor="Black",
                borderwidth=0.0,
                font=dict(size=10),
                xanchor="right",  # Anchors the legend's x position
                yanchor="top",  # Anchors the legend's y position
            ),
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
            # margin={"r": 10, "t": 40, "b": 15, "l": 50},
            margin={"r": 20, "t": 40, "b": 0, "l": 80},
            showlegend=False,
            height=300,
            width=750,
        )
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
        fig_json = fig.to_json()
        with open(json_file_path, "w") as f:
            f.write(fig_json)
        # Read the JSON file
        with open(json_file_path, "r") as json_file:
            data = json.load(json_file)
        combined_data = {"plotlyData": data, "dfData": df_sliced_json}
        fig = ""
        return combined_data

    def downloadTimesSeries(
        self,
        file_path: str,
        raster: xr.Dataset,
        downloadOption: str,
        agg_method: str,  # e.g., "mean", "max", "min"
        predicate: str,  # e.g., ">", "<", "==", "!=", ">=", "<="
        value: float,
    ):

        operator_mapping = {
            ">": operator.gt,
            "<": operator.lt,
            "=": operator.eq,  # Note that in Python, '==' is used for equality check
            "<=": operator.le,
            ">=": operator.ge,
            "!=": operator.ne,
        }
        # Aggregate over time dimension
        op_func = operator_mapping[predicate]
        variable = raster[self.variable]

        if predicate == "any":
            print(downloadOption)
            if downloadOption == "Areas":
                condition_test = op_func(variable, value)
                raster["condition_met"] = (
                    raster.where(condition_test)[self.variable]
                    .notnull()
                    .any(dim="time")
                    .astype(bool)
                )
            elif downloadOption == "Timees":
                condition_test = op_func(variable, value)
                raster["condition_met"] = (
                    raster.where(condition_test)[self.variable]
                    .notnull()
                    .any(dim=["latitude", "longitude"])
                    .astype(bool)
                )
            else:

                raise ValueError("iHARPV: Please specify data type for download")

        else:
            if downloadOption == "Areas":
                condition_test = op_func(variable, value)
                raster["condition_met"] = (
                    raster.where(condition_test)[self.variable]
                    .notnull()
                    .all(dim="time")
                    .astype(bool)
                )
            elif downloadOption == "Timees":
                condition_test = op_func(variable, value)
                raster["condition_met"] = (
                    raster.where(condition_test)[self.variable]
                    .notnull()
                    .all(dim=["latitude", "longitude"])
                    .astype(bool)
                )
            else:
                raise ValueError("iHARPV: Please specify data type for download")

        # print(raster)
        raster.to_netcdf(file_path)
        return
