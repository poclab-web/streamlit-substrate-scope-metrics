import json
import streamlit as st
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="Links")
st.title("Links")

with open("explanation/explanation.json", mode="r") as f:
    explanation = json.load(f)


st.subheader(explanation[st.session_state["language"]]["links"]["1"])
st.write("Ichizawa, K.; Nishii, T.; Gotoh, H. Substrate Scope in Organic Reactions: Comparison of Fingerprint-Based and Reactivity-Based Similarity Metrics. *Journal of Chemical Information and Modeling* **2026**, *66* (7), 3878–3891.")
st.page_link(
    "https://doi.org/10.1021/acs.jcim.5c03167",
    label="Article",
    icon=":material/link:" 
    )

st.subheader("POCLab GitHub Repository")
st.page_link(
    "https://github.com/poclab-web",
    label="POCLab (Gotoh lab in Yokohama National University)",
    icon=":material/link:"
    )

st.subheader(explanation[st.session_state["language"]]["links"]["2"])
st.page_link(
    "https://poclab-web.github.io/homepage/",
    label="POCLab",
    icon=":material/link:"
    )

st.subheader(explanation[st.session_state["language"]]["links"]["3"])
st.page_link(
    "https://www.ynu.ac.jp/english/",
    label="Yokohama National University official website",
    icon=":material/link:"
    )
