import tempfile
import typing
from pathlib import Path

import geofileops as gfo
import geopandas as gpd
import shapely

from spatial_utils.geospatial import ensure_projected_CRS, match_crs


def get_temp_files(*file_names):
    # Create a temporary folder to write data to disk since geofileops only works with on-disk
    # objects. Once it goes out of scope all the contents will be deleted from disk.
    temp_folder = tempfile.TemporaryDirectory()
    # Create a path for each of the filenames
    temp_paths = [str(Path(temp_folder.name, file_name)) for file_name in file_names]

    # The temp folder needs to be returned so it doesn't go out of scope immediately, since at that
    # point the contents of the directory are deleted.
    return [temp_folder] + temp_paths


def geofileops_buffer(input_gdf, distance, convert_to_projected_CRS=True):
    temp_folder, input_path, output_path = get_temp_files("input.gpkg", "output.gpkg")

    inital_crs = input_gdf.crs
    # If requested convert the data into a projected CRS (if not already)
    if convert_to_projected_CRS:
        input_gdf = ensure_projected_CRS(input_gdf)

    # Write the data to disk
    input_gdf.to_file(input_path)

    gfo.buffer(
        input_path=input_path,
        output_path=output_path,
        distance=distance,
        endcap_style=gfo.BufferEndCapStyle.FLAT,
        join_style=gfo.BufferJoinStyle.MITRE,
        force=True,
    )

    buffered = gpd.read_file(output_path)
    # TODO check if this is a cheap operation if the current and requested CRS match
    buffered.to_crs(inital_crs, inplace=True)
    return buffered


def geofileops_simplify(input_gdf, tolerence, convert_to_projected_CRS: bool = True):
    # Get the temporary file paths
    temp_folder, input_path, output_path = get_temp_files("input.gpkg", "output.gpkg")

    inital_crs = input_gdf.crs
    # If requested convert the data into a projected CRS (if not already)
    if convert_to_projected_CRS:
        input_gdf = ensure_projected_CRS(input_gdf)

    # Write the data to disk
    input_gdf.to_file(input_path)

    gfo.simplify(
        input_path=input_path,
        output_path=output_path,
        tolerance=tolerence,
        force=True,
    )

    # Read the result back in
    simplified = gpd.read_file(output_path)
    # TODO check if this is a cheap operation if the current and requested CRS match
    simplified.to_crs(inital_crs, inplace=True)
    return simplified


def geofileops_clip(
    input_gdf: gpd.GeoDataFrame,
    clip_geometry: typing.Union[gpd.GeoDataFrame, shapely.Geometry],
) -> gpd.GeoDataFrame:
    temp_folder, input_path, output_path, clip_path = get_temp_files(
        "input.gpkg", "output.gpkg", "clip.gpkg"
    )

    input_gdf.to_file(input_path)

    # Handle the fact that the clip geometry can be either a geodataframe or shapely geometry
    if isinstance(clip_geometry, gpd.GeoDataFrame):
        clip_geometry.to_file(clip_path)
    elif isinstance(clip_geometry, shapely.Geometry):
        gpd.GeoDataFrame(geometry=[clip_geometry]).to_file(clip_path)

    gfo.clip(
        input_path=input_path,
        output_path=output_path,
        clip_path=clip_path,
    )

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
    temp_folder, left_df_file, right_df_file, union_file = get_temp_files(
        "left.gpkg", "right.gpkg", "union.gpkg"
    )

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
    temp_folder, input_path, output_path = get_temp_files("input.gpkg", "output.gpkg")

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
