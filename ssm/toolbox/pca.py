import pandas as pd
from sklearn.decomposition import PCA


def pca_2d(df):
    pca = PCA(n_components=2, whiten=True, random_state=0)
    pca.fit(df)
    df_score = pd.DataFrame(pca.transform(df), columns=["PC1", "PC2"])
    ratio = pca.explained_variance_ratio_

    return df_score, ratio
