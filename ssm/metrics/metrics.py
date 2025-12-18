"""
Calculate metrics of the substrate scope.
"""

import pandas as pd
from scipy.spatial import ConvexHull
from scipy.spatial import distance


def metrics_range(df):
    """Calculates the sum of differences between max and min values of each descriptors.

    Returns
    - Sum of differences between max and min values of each descriptors.
    - DataFrame indicating max value, min value, and range of each descriptors.

    Args:
        df (DataFrame): DataFrame consisting only of descriptors.
    Returns:
        tuple: As mentioned above (float, DataFrame)
    """
    df.loc["min"] = df.min()
    df.loc["max"] = df.max()

    df2 = df.loc["min":"max"].transpose()
    df2["range"] = df2["max"] - df2["min"]
    range_sum = df2["range"].sum()

    return range_sum, df2


def distance_matrix(df):
    """Calculates distances between substrates and creates distance matrix.

    Args:
        df (DataFrame): DataFrame consisting only of InChI and descriptors.

    Returns:
        DataFrame: Distance matrix.
    """
    df_descriptors = df.iloc[:, 1:]
    distances = distance.pdist(df_descriptors, metric="euclidean")
    df_distance_matrix = pd.DataFrame(distance.squareform(distances))

    return pd.concat([df.iloc[:, 0], df_distance_matrix], axis=1)


def metrics_distance(df):
    """Calculates average and max values of distances from average coordinate to each substrate.

    Returns
    - Average coordinate
    - Sum of distances from average coodinate to each substrate
    - Distance_avg
    - InChI of the substrate nearest to the average coordinate
    - Distance_max
    - InChI of the substrate farthest to the average coordinate
    - Distance matrix

    Args:
        df (DataFrame): DataFrame consisting only of InChI and descriptors.

    Returns:
        tuple: As mentioned above (tuple, float, float, str, float, str, DataFrame)
    """
    name = df.columns[0]
    average = df.mean(numeric_only=True)
    average_coordinate = tuple(average.to_list())
    
    average[name] = "Average"
    average_index = len(df)
    df.loc[average_index] = average

    df2 = distance_matrix(df)
    df2["Sum"] = df2.sum(axis=1, numeric_only=True)

    distance_sum = df2.at[average_index, "Sum"]
    distance_avg = distance_sum / (len(df2) - 1)

    df3 = df2.drop(average_index)

    nearest_inchi = df3.at[df3[average_index].idxmin(), name]
    distance_max = df3[average_index].max()
    farthest_inchi = df3.at[df3[average_index].idxmax(), name]
    
    return (average_coordinate, distance_sum, distance_avg,
            nearest_inchi, distance_max, farthest_inchi, df2)


def metrics_ConvexHull(array):
    """Calculates the area and vertices of the convex hull.

    Args:
        array (array): 2-dimentional coodinates of the dataset.

    Returns:
        tuple: The area and indices of vertices.
    """
    hull = ConvexHull(array)
    area = hull.area
    vertices = hull.vertices

    return area, vertices
