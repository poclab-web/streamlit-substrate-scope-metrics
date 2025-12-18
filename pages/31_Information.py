import json
import platform
import sys

import bokeh
import numpy as np
import pandas as pd
import psutil
from rdkit import rdBase
import scipy
import sklearn
import streamlit as st

from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="Information")
st.title("Information")

with open("explanation/explanation.json", mode="r") as f:
    explanation = json.load(f)


st.header(explanation[st.session_state["language"]]["information"]["1"])

st.write(f"**{platform.platform()}**")
st.write(f"**Python**: {sys.version}")

st.space("small")

st.write(f"**Streamlit**: {st.__version__}")
st.write(f"**NumPy**: {np.__version__}")
st.write(f"**SciPy**: {scipy.__version__}")
st.write(f"**pandas**: {pd.__version__}")
st.write(f"**RDKit**: {rdBase.rdkitVersion}")
st.write(f"**scikit-learn**: {sklearn.__version__}")
st.write(f"**Bokeh**: {bokeh.__version__}")


st.space("small")

st.header("CPU")

st.write(f"{explanation[st.session_state["language"]]["information"]["2"]}: {psutil.cpu_count(logical=False)}")
st.write(f"{explanation[st.session_state["language"]]["information"]["3"]}: {psutil.cpu_count(logical=True)}")
st.write(f"{explanation[st.session_state["language"]]["information"]["4"]}: {psutil.cpu_freq().max / 1000} GHz")
st.write(f"{explanation[st.session_state["language"]]["information"]["5"]}: {psutil.cpu_percent()} %")


st.space("small")

st.header(explanation[st.session_state["language"]]["information"]["6"])
virtual_mem = psutil.virtual_memory()

st.write(f"{explanation[st.session_state["language"]]["information"]["7"]}: {virtual_mem.total / (1024 ** 3):.2f} GB")
st.write(f"{explanation[st.session_state["language"]]["information"]["8"]}: {virtual_mem.used / (1024 ** 3):.2f} GB")
st.write(f"{explanation[st.session_state["language"]]["information"]["9"]}: {virtual_mem.available / (1024 ** 3):.2f} GB")
st.write(f"{explanation[st.session_state["language"]]["information"]["10"]}: {virtual_mem.percent} %")
