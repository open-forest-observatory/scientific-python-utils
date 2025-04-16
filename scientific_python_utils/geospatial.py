import geopandas as gpd
import pyproj

from scientific_python_utils.constants import LAT_LON_CRS


def match_crs(
    target_gdf: gpd.GeoDataFrame, updateable_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Make sure the CRS of the two geodataframes match by projecting the second one

    Args:
        target_gdf (gpd.GeoDataFrame): The data from which the CRS is obtained
        updateable_gdf (gpd.GeoDataFrame): The data to be projected

    Raises:
        ValueError: If the data cannot be projected

    Returns:
        gpd.GeoDataFrame: The data contained in `updatedable_gdf` projected to the CRS of `target_gdf`
    """
    target_crs = target_gdf.crs
    updatable_crs = updateable_gdf.crs

    if target_crs is None and updatable_crs is not None:
        raise ValueError("Target CRS is None while the updateable CRS is not")

    # If the CRS don't match, transform the updatable one
    if target_crs != updatable_crs:
        updateable_gdf = updateable_gdf.to_crs(target_crs)

    # TODO think about how to handle if both are None - do nothing?
    # TODO consider whether a copy should be returned in the case where no update is made
    return updateable_gdf


def ensure_projected_CRS(geodata: gpd.GeoDataFrame):
    """Returns a projected geodataframe from the provided geodataframe by converting it to
    ESPG:4326 (if not already) and determining the projected CRS from the point
    coordinates.

    Args:
        geodata (gpd.GeoDataGrame): Original geodataframe that is potentially unprojected
    Returns:
        gpd.GeoDataGrame: projected geodataframe
    """
    # If CRS is projected return immediately
    if geodata.crs.is_projected:
        return geodata

    # If CRS is geographic and not long-lat, convert it to long-lat
    if geodata.crs.is_geographic and geodata.crs != LAT_LON_CRS:
        geodata = geodata.to_crs(LAT_LON_CRS)

    # Convert geographic long-lat CRS to projected CRS
    point = geodata["geometry"].iloc[0].centroid
    geometric_crs = get_projected_CRS(lon=point.x, lat=point.y)
    return geodata.to_crs(geometric_crs)


def get_projected_CRS(lat, lon, assume_western_hem=True):
    if assume_western_hem and lon > 0:
        lon = -lon
    # https://gis.stackexchange.com/questions/190198/how-to-get-appropriate-crs-for-a-position-specified-in-lat-lon-coordinates
    epgs_code = 32700 - round((45 + lat) / 90) * 100 + round((183 + lon) / 6)
    crs = pyproj.CRS.from_epsg(epgs_code)
    return crs
