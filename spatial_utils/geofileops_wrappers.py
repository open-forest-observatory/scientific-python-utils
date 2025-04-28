import tempfile
import typing
from pathlib import Path

import geofileops as gfo
import geopandas as gpd

from spatial_utils.geospatial import match_crs


def geofileops_clip(
    input_gdf: gpd.GeoDataFrame, clip_geometry: typing.Union[gpd.GeoDataFrame]
) -> gpd.GeoDataFrame:
    # Create a temporary folder to write data to disk since geofileops only works with on-disk
    # objects. Once it goes out of scope all the contents will be deleted from disk.
    temp_folder = tempfile.TemporaryDirectory()
    input_path = str(Path(temp_folder.name, "input.gpkg"))
    output_path = str(Path(temp_folder.name, "output.gpkg"))
    clip_path = str(Path(temp_folder.name, "clip.gpkg"))

    input_gdf.to_file(input_path)
    # TODO better handle other datatypes
    clip_geometry.to_file(clip_path)

    gfo.clip(input_path=input_path, output_path=output_path, clip_path=clip_path)

    clipped = gpd.read_file(output_path)
    return clipped


def geofileops_overlay(
    left_gdf: gpd.GeoDataFrame,
    right_gdf: gpd.GeoDataFrame,
    input1_columns_prefix: str = "l1_",
    input2_columns_prefix: str = "l2_",
) -> gpd.GeoDataFrame:
    """Perform an operation similar to geopandas.overlay(how="union")

    Args:
        left_gdf (gpd.GeoDataFrame):
            First object to union. Note that the right_gdf will be transformed to this CRS if the
            two CRS do not match.
        right_gdf (gpd.GeoDataFrame): Second object to union
        input1_columns_prefix (str, optional):
            Prefix to add to columns in the last dataframe. Note that this is applied to all columns,
            not just ones that colide between the two dataframes, which differs from native
            geopandas. Defaults to "l1_".
        input2_columns_prefix (str, optional): See `input1_columns_prefix`. Defaults to "l2_".

    Returns:
        gpd.GeoDataFrame: The overlaid data in the CRS of the first dataframe.
    """
    # Create a temporary folder to write data to disk since geofileops only works with on-disk
    # objects. Once it goes out of scope all the contents will be deleted from disk.
    temp_folder = tempfile.TemporaryDirectory()
    left_df_file = str(Path(temp_folder.name, "left.gpkg"))
    right_df_file = str(Path(temp_folder.name, "right.gpkg"))
    union_file = str(Path(temp_folder.name, "union.gpkg"))

    # Make the CRS match between the two datasets or error if impossible
    right_gdf = match_crs(target_gdf=left_gdf, updateable_gdf=right_gdf)

    # Write data to disk
    left_gdf.to_file(left_df_file)
    right_gdf.to_file(right_df_file)

    # This performs the same operation as geopandas overlay(how="union")
    gfo.union(
        input1_path=left_df_file,
        input2_path=right_df_file,
        output_path=union_file,
        input1_columns_prefix=input1_columns_prefix,
        input2_columns_prefix=input2_columns_prefix,
    )
    # Read the result back in
    overlay_data = gpd.read_file(union_file)

    return overlay_data


def geofileops_dissolve(
    input_gdf: gpd.GeoDataFrame,
    groupby_columns: typing.Optional[typing.Union[str, typing.List[str]]] = None,
    retain_all_columns: bool = True,
) -> gpd.GeoDataFrame:
    """Perform an operation similar to a geopandas.dissolve

    Args:
        input_gdf (gpd.GeoDataFrame): Input data
        groupby_columns (typing.Optional[typing.Union[str, typing.List[str]]], optional):
            Column or columns to dissolve by. If not specified, all data will be merged. Defaults to None.
        retain_all_columns (bool, optional):
            Should the data from columns which are not dissolved on be retained. Defaults to True.

    Returns:
        gpd.GeoDataFrame: Dissolved dataframe
    """
    # Create a temporary folder to write data to disk since geofileops only works with on-disk
    # objects. Once it goes out of scope all the contents will be deleted from disk.
    temp_folder = tempfile.TemporaryDirectory()
    input_path = Path(temp_folder.name, "input.gpkg")
    output_path = Path(temp_folder.name, "output.gpkg")

    # Write the input data to disk
    input_gdf.to_file(input_path)

    if retain_all_columns:
        # Similar to a geopandas dissolve, the minimum value across all aggregated rows is kept after
        # dissolving

        # Default, no columns provided
        if groupby_columns is None:
            groupby_columns_list = []
        # Only one column is provided
        elif not isinstance(groupby_columns, (list, tuple)):
            groupby_columns_list = [groupby_columns]
        else:
            groupby_columns_list = groupby_columns

        agg_columns = list(
            set(input_gdf.columns) - set(groupby_columns_list + ["geometry"])
        )
        agg_columns = {
            "columns": [{"column": c, "as": c, "agg": "min"} for c in agg_columns]
        }
    else:
        # All columns except for the ones that are being aggregated by will be dropped
        agg_columns = None

    # Dissolve
    gfo.dissolve(
        input_path=input_path,
        output_path=output_path,
        explodecollections=False,
        groupby_columns=groupby_columns,
        agg_columns=agg_columns,
    )

    # Read the result and return
    dissolved = gpd.read_file(output_path)
    return dissolved
