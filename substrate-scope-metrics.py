import io
import json

import pandas as pd
from sklearn.preprocessing import StandardScaler
import streamlit as st

from ssm import smiles_to_inchi
from utils.sidebar import sidebar
from utils.initialize_session_state import reset


sidebar()
st.set_page_config(page_title="Substrate Scope Metrics")
st.title("Substrate Scope Metrics")

with open("explanation/explanation.json", mode="r") as f:
    explanation = json.load(f)

st.write(explanation[st.session_state["language"]]["home"]["1"])
st.space(size="small")


how_to_use = st.expander(explanation[st.session_state["language"]]["home"]["2"])
with open(f"explanation/how_to_use_{st.session_state['language']}.md", mode="r") as f:
    how_to_use.write(f.read())

st.space(size="small")


if st.session_state["loaded"] == False:
    st.header(explanation[st.session_state["language"]]["home"]["3"])

    if st.toggle(explanation[st.session_state["language"]]["home"]["4"], value=True):
        st.session_state["standardize"] = True
    else:
        st.session_state["standardize"] = False

    st.session_state["n_conditions"] = st.number_input(
        explanation[st.session_state["language"]]["home"]["5"], min_value=1, value=None
        )

    dataset = st.file_uploader(
        explanation[st.session_state["language"]]["home"]["6"], type="csv"
        )

    if dataset is not None:
        if st.session_state["n_conditions"] is None:
            st.warning(
                explanation[st.session_state["language"]]["home"]["7"], icon=":material/warning:"
                )
        else:
            with st.spinner(explanation[st.session_state["language"]]["common"]["5"]):
                st.session_state["dataset"] = dataset.getvalue().decode("utf-8")
                st.session_state["dataset_name"] = dataset.name

                df1 = pd.read_csv(io.StringIO(st.session_state["dataset"]))
                df1.iloc[:, 0] = smiles_to_inchi(df1.iloc[:, 0])
                st.session_state["dataset"] = df1.to_csv(index=False)

                st.session_state["loaded"] = True
                st.rerun()

else:
    st.header(explanation[st.session_state["language"]]["home"]["8"])


if st.session_state["dataset"] is not None:
    df1 = pd.read_csv(io.StringIO(st.session_state["dataset"]))

    if st.session_state["standardize"]:
        if st.session_state["standardized"] == False:
            df2 = df1.iloc[:, (st.session_state["n_conditions"] + 1):]
            
            if df2.shape[1] == 0:
                st.session_state["standardize"] = False
                st.rerun()

            standardized_df2 = StandardScaler().fit_transform(df2)
            df_standardized = pd.DataFrame(standardized_df2, columns=df2.columns)

            df3 = pd.concat([df1.iloc[:, :(st.session_state["n_conditions"] + 1)],
                             df_standardized],
                             axis=1)

            st.session_state["dataset"] = df3.to_csv(index=False)
            st.session_state["standardized"] = True

    with st.container(border=True):
        st.subheader(st.session_state["dataset_name"])
        st.write(
            f"**{explanation[st.session_state['language']]['home']['9']}**: {len(df1)}"
            )
        st.write(
            f"**{explanation[st.session_state['language']]['home']['10']}**: {len(df1.columns) - st.session_state['n_conditions'] - 1}"
            )

        if st.session_state["standardized"]:
            st.badge(
                explanation[st.session_state["language"]]["home"]["11"],
                icon=":material/check:",
                color="primary"
                )
            
        with st.expander(explanation[st.session_state["language"]]["home"]["17"]):
            for descriptor in df1.columns[st.session_state["n_conditions"] + 1:]:
                st.write(descriptor)
        
        st.divider()

        for i in range(st.session_state["n_conditions"]):
            condition = df1.columns[i+1]
            st.subheader(condition)

            count = len(df1[df1[condition] == 1])
            st.write(
                f"**{explanation[st.session_state['language']]['common']['1']}**: {count}"
                )

            n_negative = len(df1[df1[condition] == -1])
            if n_negative > 0:
                st.write(
                    f"**{explanation[st.session_state['language']]['common']['2']}**: {n_negative}"
                    )

            if len(st.session_state["conditions"]) != st.session_state["n_conditions"]:
                st.session_state["conditions"].append(condition)
                st.session_state["label"].append(f"{condition} ({count})")


    with st.expander(
        explanation[st.session_state["language"]]["home"]["12"],
        icon=":material/settings:"
        ):
        st.subheader("ECFP")
        ecfp_settings_0 = st.number_input(
            explanation[st.session_state["language"]]["home"]["13"],
            value=st.session_state["ecfp_settings"][0]
            )
        ecfp_settings_1 = st.number_input(
            explanation[st.session_state["language"]]["home"]["14"],
            value=st.session_state["ecfp_settings"][1]
            )
        
        if st.button(
            explanation[st.session_state["language"]]["home"]["15"],
            width="stretch"
            ):
            st.session_state["ecfp_settings"] = [ecfp_settings_0, ecfp_settings_1]
            st.rerun()

        
    st.button(
        explanation[st.session_state["language"]]["home"]["16"],
        on_click=reset,
        width="stretch"
        )
