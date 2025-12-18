import io
import json

import numpy as np
import pandas as pd
import streamlit as st

from ssm import metrics_range
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="Range")
st.title("Range")

with open("explanation/explanation.json", mode="r") as f:
    explanation = json.load(f)


if st.session_state["dataset"] is None:
    st.warning(
        explanation[st.session_state["language"]]["common"]["3"],
        icon=":material/warning:"
        )

else:
    with st.spinner(explanation[st.session_state["language"]]["common"]["4"]):
        df1 = pd.read_csv(io.StringIO(st.session_state["dataset"]))

        li_range = []
        li_details = []
        for condition in st.session_state["conditions"]:
            df2 = df1[df1[condition] == 1].iloc[:, (st.session_state["n_conditions"] + 1):].copy()

            ss = metrics_range(df2)
            li_range.append(ss[0])
            li_details.append(ss[1]["range"])

        ss = metrics_range(df1.iloc[:, (st.session_state["n_conditions"] + 1):].copy())
        li_range.append(ss[0])

        df_details = pd.concat(li_details, axis=1)
        df_details = df_details.set_axis(
            st.session_state["conditions"], axis=1, copy=False
            )


    df_result = pd.DataFrame(
        {"Condition": st.session_state["label"] + [f"Dataset ({len(df1)})"],
         "Range": li_range}
        )
    
    df_coverage = pd.DataFrame(
        {"Condition": st.session_state["label"],
         "Range (Coverage)": np.array(li_range[:-1]) / li_range[-1] * 100}
        )
    
    
    tab1, tab2 = st.tabs([
        explanation[st.session_state["language"]]["common"]["7"],
        explanation[st.session_state["language"]]["common"]["8"]
    ])

    with tab1:
        st.bar_chart(
            df_result.iloc[:-1],
            x="Condition",
            y="Range",
            y_label="",
            horizontal=True,
            sort=False
            )
        
    with tab2:
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="Range (Coverage)",
            x_label="Range (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )

    st.divider()


    st.header(explanation[st.session_state["language"]]["common"]["6"])
    st.subheader(explanation[st.session_state["language"]]["range"]["3"])

    color = None
    d28 = [
        "m_1_L", "m_1_B1", "m_1_B5",
        "p_L", "p_B1", "p_B5",
        "m_2_L", "m_2_B1", "m_2_B5",
        "o_L", "o_B1", "o_B5",
        "Charge_O", "Charge_next_O", "Charge_o_1", "Charge_m_1", "Charge_p", "Charge_m_2", "Charge_o_2",
        "HOMO", "MolLogP", "MolWt",
        "num_N", "num_O", "num_S", "num_F", "num_Cl", "num_Br"
        ]
    
    if list(df1.iloc[:, (st.session_state["n_conditions"] + 1):].columns) == d28:
        color = [
            "#FBEED0", "#F2CC76", "#F6AA00",
            "#D3EEE4", "#82CCB1", "#03AF7A",
            "#E5CEE8", "#B56DBC", "#990099",
            "#DFF2FE", "#A3DAFB", "#4DC4FF",
            "#FCEDE6", "#F8D2C2", "#F5BA9E", "#F2A07B", "#EF8658", "#EF6F3B", "#FF4B00",
            "#CFDEFC", "#709DF8", "#005AFF",
            "#FBE6E6", "#F4B5B5", "#FF8082", "#E3D8CD", "#AB8D6B", "#804000"
            ]


    tab1, tab2 = st.tabs([
        explanation[st.session_state["language"]]["range"]["1"],
        explanation[st.session_state["language"]]["range"]["2"]
    ])

    with tab1:
        st.bar_chart(
            df_details.T,
            color=color,
            horizontal=True,
            sort=False,
            height=500
            )
        
    with tab2:
        st.dataframe(df_details)
    

    st.divider()
    
    st.subheader("Range")
    st.dataframe(df_result, hide_index=True)
    