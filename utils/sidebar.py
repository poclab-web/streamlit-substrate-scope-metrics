import streamlit as st
from utils.initialize_session_state import initialize_session_state


def sidebar():
    initialize_session_state()

    st.logo("logo/poclab-logo-white-1.png", icon_image="logo/poclab-logo-white-2.png")

    languages = {
        "English": "en",
        "日本語": "ja"
    }

    with st.sidebar:
        st.write("Substrate Scope Metrics")

        st.page_link("substrate-scope-metrics.py", label="Home", icon=":material/home:")

        st.space("small")
        
        st.page_link("pages/11_ECFP_OnBits.py", label="ECFP OnBits", icon=":material/switches:")
        st.page_link("pages/12_ECFP_Tanimoto.py", label="ECFP Tanimoto", icon=":material/fingerprint:")
        st.page_link("pages/13_ECFP_ConvexHull.py", label="ECFP ConvexHull", icon=":material/bubble_chart:")

        st.space("small")
        
        st.page_link("pages/21_Range.py", label="Range", icon=":material/arrow_range:")
        st.page_link("pages/22_Distance.py", label="Distance", icon=":material/arrows_output:")
        st.page_link("pages/23_ConvexHull.py", label="ConvexHull", icon=":material/scatter_plot:")

        st.space("small")
        
        st.page_link("pages/31_Information.py", label="Information", icon=":material/info:")
        st.page_link("pages/32_Links.py", label="Links", icon=":material/link:")

        st.divider()
        
        select_language = st.selectbox(
            ":material/language: Language",
            languages.keys(),
            index=list(languages.values()).index(st.session_state["language"])
            )
        if st.session_state["language"] != languages[select_language]:
            st.session_state["language"] = languages[select_language]
            st.rerun()
