import io
from multiprocessing import Pool, cpu_count

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import streamlit as st


def ecfp(inchi_list, radius, fpSize):
    """Generates ECFP.

    Args:
        inchi_list (array-like): List of InChI.
        radius (int, optional): The number of iterations to grow the fingerprint.
        fpSize (int, optional): Size of the generated fingerprint.
    Returns:
        tuple: A tuple of ExplicitBitVects.
    """
    with Pool(cpu_count()) as p:
        molecules = p.map(Chem.MolFromInchi, inchi_list)

    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize)
    fps = mfpgen.GetFingerprints(molecules, numThreads=cpu_count())

    return fps


def generate_ecfp():
    if st.session_state["generated_ecfp"] != st.session_state["ecfp_settings"]:
        df = pd.read_csv(io.StringIO(st.session_state["dataset"]))

        st.session_state["ecfp"] = ecfp(
            df.iloc[:, 0],
            st.session_state["ecfp_settings"][0],
            st.session_state["ecfp_settings"][1]
            )
        
        st.session_state["generated_ecfp"] = st.session_state["ecfp_settings"]
        st.session_state["ecfp_pca"] = None

    else:
        pass
