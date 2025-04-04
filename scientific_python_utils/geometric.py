import numpy as np
import geopandas as gpd

def merge_classified_polygons(classified_polygons, class_column=None):
    grouped = {k[0]:v for k, v in classified_polygons.groupby([class_column])}

    unique_classes = list(grouped.keys())

    counts_overlays = []

    for k, v in grouped.items():
        v.drop("classes", axis=1, inplace=True)
        v["counts"] = 1

        overlay = v.iloc[0:1]

        for i in range(1, len(v)):
            overlay = gpd.overlay(overlay, v.iloc[i:i+1], how="union")
            overlay.fillna(0, inplace=True)
            overlay["counts"] = overlay["counts_1"] + overlay["counts_2"]
            overlay.drop(["counts_1", "counts_2"], axis=1, inplace=True)
        overlay = overlay.dissolve("counts", as_index=False)
        overlay[k] = overlay["counts"]
        overlay.drop("counts", axis=1, inplace=True)

        counts_overlays.append(overlay)

    votes_per_class = counts_overlays[0]
    for single_class_overlay in counts_overlays[1:]:
        votes_per_class = gpd.overlay(votes_per_class, single_class_overlay, how="union")

    votes_per_class.fillna(0, inplace=True)

    class_counts_matrix = votes_per_class[unique_classes].values
    max_class_counts = np.max(class_counts_matrix, axis=1, keepdims=True)

    # Find rows where one class has the most votes
    one_max_class = np.sum(max_class_counts == class_counts_matrix, axis=1) == 1

    # Extract rows where one class has the most votes
    rows_with_one_class = votes_per_class.iloc[one_max_class]
    # Label them with the max class
    rows_with_one_class["max_class"] = rows_with_one_class[unique_classes].idxmax(axis=1)
    # Dissolve all polygons so we have one (multi)polygon per class
    rows_with_one_class = rows_with_one_class.dissolve("max_class", as_index=False)
    rows_with_one_class.plot("max_class", cmap="tab10")

    # Compute the area of each
    rows_with_one_class["area"] = rows_with_one_class.area

    sorted_inds = (rows_with_one_class["area"]).argsort()
    sorted_classes = rows_with_one_class.index[sorted_inds].tolist()

    # Determine which classes (if any) have no non-overlapping regions. Add them to the start of the
    # list
    zero_area_classes = [c for c in unique_classes if c not in sorted_classes]
    # Prepend the classes to the beginning of the list of sorted classes
    sorted_classes = zero_area_classes + sorted_classes

    # Get the new order of column names
    max_class = votes_per_class[sorted_classes].idxmax(axis=1)
    votes_per_class["max_class"] = max_class
    votes_per_class.plot("max_class")