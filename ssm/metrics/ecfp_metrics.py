from itertools import combinations
from multiprocessing import Pool, cpu_count

import numpy as np
from rdkit.Chem import DataStructs


def fp_dropna(df_ecfp):
    """Deletes the bits taking 0 in all molecules.

    Args:
        df_ecfp (DataFrame): DataFrame of ECFP.

    Returns:
        DataFrame: DataFrame of ECFP without the bit taking 0 in all molecules.
    """
    df_ecfp.loc[len(df_ecfp)] = df_ecfp.sum(numeric_only=True)
    df2 = df_ecfp.T[df_ecfp.T[len(df_ecfp) - 1] != 0].T
    
    return df2


def on_bits(df_ecfp):
    """Calculates ECFP OnBits.

    Args:
        df_ecfp (DataFrame): DataFrame of ECFP.

    Returns:
        tuple: ECFP OnBits (int), Bit numbers (index), DataFrame
    """
    df2 = fp_dropna(df_ecfp)
    bit_list = df2.columns
    ans = len(bit_list)

    return ans, bit_list, df2


def tanimoto_distance(pair):
    return 1 - DataStructs.TanimotoSimilarity(pair[0], pair[1])


def tanimoto(df_ecfp_vect, get_indices=True):
    """Calculates Tanimoto distances between molecules.

    Tanimoto distance = 1 - Tanimoto similarity

    Returns
    - Average value
    - Maximum value
    - Minimum value
    - Combination of molecule indices providing maximum value
    - Combination of molecule indices providing minimum value

    Args:
        df_ecfp (DataFrame): DataFrame of ECFP as ExplicitBitVects.
        get_indices (bool, optional): Whether to return combination of molecule indices providing max or min value. Defaults to True.

    Returns:
        tuple: As mentioned above (np.float64, float, float, list, list)
    """
    comb = combinations(df_ecfp_vect.iloc[:, 0], 2)

    if len(df_ecfp_vect) > 2000:
        with Pool(cpu_count()) as p:
            distances = p.map(tanimoto_distance, comb)
    
    else:
        distances = list(map(tanimoto_distance, comb))

    distance_avg = np.mean(distances)
    distance_max = max(distances)
    distance_min = min(distances)

    max_indices = []
    min_indices = []
    if get_indices:
        comb_index = list(combinations(df_ecfp_vect.index, 2))
        for i, v in enumerate(distances):
            if v == distance_max:
                max_indices.append([comb_index[i][0], comb_index[i][1]])
            elif v == distance_min:
                min_indices.append([comb_index[i][0], comb_index[i][1]])

    return distance_avg, distance_max, distance_min, max_indices, min_indices
