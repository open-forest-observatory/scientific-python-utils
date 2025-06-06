from typing import List, Optional, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from geopandas import GeoDataFrame, GeoSeries
from shapely import (
    Geometry,
    MultiPolygon,
    Polygon,
    difference,
    intersection,
    make_valid,
    union,
)
from tqdm import tqdm

from spatial_utils.geofileops_wrappers import geofileops_dissolve, geofileops_overlay


def merge_classified_polygons_by_voting(
    classified_polygons: gpd.GeoDataFrame,
    class_column: str,
    tiebreaking_class_order: Optional[List[str]] = None,
    print_tiebreaking_stats: bool = False,
) -> gpd.GeoDataFrame:
    """
    Take multiple potentially-overlapping polygons with associated class information and merge them
    into non-overlapping classified polygons using the following strategy.
    * For each location, compute the number of polygons that cast votes for each class
    * For any regions that have a non-tied number of votes, set the output class to the majority vote
    * For any regions with tied votes, break ties in the order specified in `tiebreaking_class_order`
      if specified, or in favor of the less common class (based on un-ambigous regions)

    Args:
        classified_polygons (gpd.GeoDataFram): A geodataframe containing the `class_column` column
        class_column (str): The column to use as the class
        tiebreaking_class_order (List[str], optional):
            A list of ordered class names. Ties will be broken in favor of classes earlier in the list.
        print_tiebreaking_stats (bool, optional): Print the fraction of area that needed to be tiebroken

    Returns:
        gpd.GeoDataFrame:
            A dataframe with only the `class_column` column representing the merged class

    """
    # Create a dictionary where the keys are the classes and the values are the dataframe with all
    # the (multi)polygons corresponding to that class
    grouped_by_class = {k[0]: v for k, v in classified_polygons.groupby([class_column])}
    unique_classes = list(grouped_by_class.keys())
    # Now, for each individual class, compute how many overlapping polygons there are for each
    # location. This information will be used in future steps to prioritize that class if there are
    # multiple predictions for the same location.
    n_overlapping_per_class = []
    for cls, geometries_per_class in grouped_by_class.items():
        # Create a new column for counts that is all ones
        geometries_per_class["counts"] = 1
        # Only retain the geometry column and the counts
        geometries_per_class = geometries_per_class[["geometry", "counts"]]

        # Initialize the number overlapping to the first geometry, which has a count of 1
        n_overlapping = geometries_per_class.iloc[0:1]
        # Iterate over the subsequent geometries
        for i in range(1, len(geometries_per_class)):
            geometry = geometries_per_class[i : i + 1]
            # Overlay the running total of overlaps with the new geometry
            n_overlapping = geofileops_overlay(n_overlapping, geometry)
            # Since we use the union approach, there will be some nan values for counts where one
            # side of the overlay did not overlap the other one. Set these to zero.
            n_overlapping.fillna(0, inplace=True)
            # Add the counts from the left and right dataframe and set them to a new "counts" column
            n_overlapping["counts"] = (
                n_overlapping["l1_counts"] + n_overlapping["l2_counts"]
            )
            # Since we'll be repeating this process, it's important to drop these columns or there
            # will be a name collision on the next iteration
            n_overlapping.drop(["l1_counts", "l2_counts"], axis=1, inplace=True)

        # TODO, consider whether it would be more efficient to dissolve after every operation. It
        # would be an additional step, but it could speed up the overlay
        n_overlapping = geofileops_dissolve(n_overlapping, groupby_columns="counts")
        # Rename the counts column to the name of the class
        n_overlapping.rename(columns={"counts": str(cls)}, inplace=True)
        # Append to the running list
        n_overlapping_per_class.append(n_overlapping)
    print("Done with the n-overlapping steps")
    # Overlay the votes for each class to get all the regions with distinct combinations of votes
    votes_per_class = n_overlapping_per_class[0]
    for i, single_class_overlay in enumerate(n_overlapping_per_class[1:]):
        votes_per_class = geofileops_overlay(
            votes_per_class,
            single_class_overlay,
            input1_columns_prefix="",
            input2_columns_prefix="",
        )
    # Similar to before, since we're doing a "union" overlay, there will be rows that don't have
    # values for all columns. Fill them in with zero.
    votes_per_class.fillna(0, inplace=True)
    # Extract the counts columns to a numpy array and find the most common class per row
    unique_classes_str = [str(x) for x in unique_classes]
    class_counts_matrix = votes_per_class[unique_classes_str].values
    max_class_counts = np.max(class_counts_matrix, axis=1, keepdims=True)

    # Find rows where one class has the most votes (there are no ties)
    one_max_class = np.sum(max_class_counts == class_counts_matrix, axis=1) == 1

    # Extract rows where one class has the most votes
    rows_with_one_class = votes_per_class.iloc[one_max_class]
    # Label them with the max class
    rows_with_one_class["max_class"] = rows_with_one_class[unique_classes_str].idxmax(
        axis=1
    )
    # Dissolve all polygons so we have one (multi)polygon per class
    rows_with_one_class = geofileops_dissolve(
        rows_with_one_class, groupby_columns="max_class"
    )

    # Compute the area of each
    rows_with_one_class["area"] = rows_with_one_class.area

    if tiebreaking_class_order is None:
        # Order the classes from smallest area to largest, based on unambigous regions
        sorted_inds = (rows_with_one_class["area"]).argsort()
        tiebreaking_class_order = rows_with_one_class["max_class"][sorted_inds].tolist()

        # Determine which classes (if any) have no non-overlapping regions. Add them to the start of the
        # list
        zero_area_classes = [
            c for c in unique_classes_str if c not in tiebreaking_class_order
        ]
        # Prepend the classes to the beginning of the list of sorted classes
        tiebreaking_class_order = zero_area_classes + tiebreaking_class_order
    else:
        breakpoint()
        tiebreaking_class_order = list(
            filter(
                lambda x: x in votes_per_class.columns.to_list(),
                tiebreaking_class_order,
            )
        )

    if print_tiebreaking_stats:
        area_of_sorted = geofileops_dissolve(rows_with_one_class).area[0]
        total_area = geofileops_dissolve(votes_per_class).area[0]

        print(
            f"Ties had to be broken for {(100 *(1 - (area_of_sorted/total_area))):.1f}% of the total predictions"
        )

    # Reorder the columns starting with the rarest class. Then compute the index of the max value.
    # Since this returns the first instance of the maximum value, ties will be broken in favor of
    # the class that had the least area in the unambigious region.
    max_class = votes_per_class[tiebreaking_class_order].idxmax(axis=1)
    # Create a new column for the max class
    votes_per_class[class_column] = max_class
    votes_per_class = votes_per_class[[class_column, "geometry"]]
    # Dissolve so there's only one (multi)polygon per class
    votes_per_class = geofileops_dissolve(votes_per_class, groupby_columns=class_column)

    # Convert the class_column values back from the str representation to the original type
    max_classes_str = votes_per_class[class_column].to_list()
    # Get the index in the list of unique values
    max_classes_index_in_unique = [unique_classes_str.index(c) for c in max_classes_str]
    # The unique values are ordered the same for the string-typed and original type. Extract the
    # corresponding value in the original type.
    max_classes_original_type = [unique_classes[i] for i in max_classes_index_in_unique]
    # Update the column
    votes_per_class[class_column] = max_classes_original_type

    return votes_per_class


def ensure_non_overlapping_polygons(
    geometries: Union[List[Geometry], gpd.GeoDataFrame],
    inplace: bool = False,
):
    # Make sure geometries is a list of shapely objects
    if isinstance(geometries, gpd.GeoDataFrame):
        original_gdf = geometries
        geometries = geometries.geometry.tolist()
    else:
        original_gdf = None

    output_geometries = [None] * len(geometries)
    union_of_added_geoms = MultiPolygon()

    areas = [geom.area for geom in geometries]
    sorted_inds = np.argsort(areas)

    for ind in tqdm(sorted_inds):
        # Get the input geometry and ensure it's valid
        geom = make_valid(geometries[ind])
        # Subtract the union of all
        geom_to_add = difference(geom, union_of_added_geoms)
        output_geometries[ind] = geom_to_add
        # Add the original geom, not the difference'd one, to avoid boundary artifacts
        union_of_added_geoms = union(geom, union_of_added_geoms)

    if original_gdf is None:
        return output_geometries
    elif inplace:
        original_gdf.geometry = output_geometries
    else:
        output_gdf = original_gdf.copy()
        output_gdf.geometry = output_geometries
        return output_gdf


def find_union_of_intersections(list_of_multipolygons, crs, vis=False):
    all_intersections = MultiPolygon()
    for i, multipolygon_a in enumerate(list_of_multipolygons):
        for multipolygon_b in list_of_multipolygons[:i]:
            new_intersection = intersection(multipolygon_a, multipolygon_b)
            all_intersections = union(all_intersections, new_intersection)
    if vis:
        geopandas_all_intersections = GeoDataFrame(
            geometry=[all_intersections], crs=crs
        )
        geopandas_all_intersections.plot()
        plt.show()
    return all_intersections


def intersects_union_of_polygons(
    query_polygons: GeoDataFrame,
    region_polygon: Union[GeoDataFrame, GeoSeries, Polygon, MultiPolygon],
):
    if isinstance(region_polygon, GeoDataFrame):
        region_polygon.plot()
        # Try to make geometries valid
        region_polygon.geometry = region_polygon.buffer(0)
        region_polygon = region_polygon.dissolve()
        region_polygon = region_polygon.geometry[0]

    # Find the polygons that are within the bounds of the raster
    intersection = query_polygons.intersection(region_polygon)
    empty_geometry = intersection.is_empty.to_numpy()
    within_bounds_IDs = np.where(np.logical_not(empty_geometry))[0]
    return within_bounds_IDs
