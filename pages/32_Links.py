import json
import streamlit as st
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="Links")
st.title("Links")

with open("explanation/explanation.json", mode="r") as f:
    explanation = json.load(f)


st.subheader(explanation[st.session_state["language"]]["links"]["1"])
st.page_link(
    "https://doi.org/",
    label="in preparation"
    )

st.subheader("POCLab GitHub Repository")
st.page_link(
    "https://github.com/poclab-web",
    label="POCLab (Gotoh lab in Yokohama National University)"
    )

st.subheader(explanation[st.session_state["language"]]["links"]["2"])
st.page_link(
    "https://poclab-web.github.io/homepage/",
    label="POCLab"
    )

st.subheader(explanation[st.session_state["language"]]["links"]["3"])
st.page_link(
    "https://www.ynu.ac.jp",
    label="Yokohama National University official website"
    )
