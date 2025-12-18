import streamlit as st


def initialize_session_state():

    # sidebar
    if "language" not in st.session_state:
        st.session_state["language"] = "en"
    
    # ecfp
    if "ecfp" not in st.session_state:
        st.session_state["ecfp"] = None

    if "generated_ecfp" not in st.session_state:
        st.session_state["generated_ecfp"] = [None, None]

    # pca
    if "pca" not in st.session_state:
        st.session_state["pca"] = None

    if "pca_ratio" not in st.session_state:
        st.session_state["pca_ratio"] = None

    if "ecfp_pca" not in st.session_state:
        st.session_state["ecfp_pca"] = None

    if "ecfp_pca_ratio" not in st.session_state:
        st.session_state["ecfp_pca_ratio"] = None

    # home
    if "loaded" not in st.session_state:
        st.session_state["loaded"] = False

    if "standardize" not in st.session_state:
        st.session_state["standardize"] = True

    if "n_conditions" not in st.session_state:
        st.session_state["n_conditions"] = None

    if "dataset" not in st.session_state:
        st.session_state["dataset"] = None

    if "dataset_name" not in st.session_state:
        st.session_state["dataset_name"] = None

    if "standardized" not in st.session_state:
        st.session_state["standardized"] = False

    if "conditions" not in st.session_state:
        st.session_state["conditions"] = []

    if "label" not in st.session_state:
        st.session_state["label"] = []

    if "ecfp_settings" not in st.session_state:
        st.session_state["ecfp_settings"] = [2, 2048]


def reset():
    for key in st.session_state.keys():
        del st.session_state[key]
    
    st.cache_data.clear()
