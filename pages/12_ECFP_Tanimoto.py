import io
import json

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st

from ssm import generate_ecfp, tanimoto
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="ECFP Tanimoto")
st.title("ECFP Tanimoto")

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

        generate_ecfp()
        df_ecfp_vect = pd.DataFrame(st.session_state["ecfp"])

        li_tanimoto_avg = []
        li_tanimoto_max = []
        li_tanimoto_min = []
        li_max_indices = []
        li_min_indices = []
        for condition in st.session_state["conditions"]:
            df2_ecfp_vect = df_ecfp_vect[df1[condition] == 1].copy()

            ss = tanimoto(df2_ecfp_vect)
            li_tanimoto_avg.append(ss[0])
            li_tanimoto_max.append(ss[1])
            li_tanimoto_min.append(ss[2])
            li_max_indices.append(ss[3])
            li_min_indices.append(ss[4])

        ss = tanimoto(df_ecfp_vect, get_indices=False)
        li_tanimoto_avg.append(ss[0])
        li_tanimoto_max.append(ss[1])
        li_tanimoto_min.append(ss[2])


    df_result = pd.DataFrame(
        {"Condition": st.session_state["label"] + [f"Dataset ({len(df1)})"],
         "ECFP Tanimoto_avg": li_tanimoto_avg,
         "ECFP Tanimoto_max": li_tanimoto_max,
         "ECFP Tanimoto_min": li_tanimoto_min}
        )
    
    df_coverage = pd.DataFrame(
        {"Condition": st.session_state["label"],
         "ECFP Tanimoto_avg (Coverage)": np.array(li_tanimoto_avg[:-1]) / li_tanimoto_avg[-1] * 100,
         "ECFP Tanimoto_max (Coverage)": np.array(li_tanimoto_max[:-1]) / li_tanimoto_max[-1] * 100}
        )
    

    tab1, tab2 = st.tabs([
        explanation[st.session_state["language"]]["common"]["7"],
        explanation[st.session_state["language"]]["common"]["8"]
    ])

    with tab1:
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["1"])
        st.bar_chart(
            df_result.iloc[:-1],
            x="Condition",
            y="ECFP Tanimoto_avg",
            y_label="",
            horizontal=True,
            sort=False
            )
        
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["2"])
        st.bar_chart(
            df_result.iloc[:-1],
            x="Condition",
            y="ECFP Tanimoto_max",
            y_label="",
            horizontal=True,
            sort=False
            )
        
    with tab2:
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["1"])
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="ECFP Tanimoto_avg (Coverage)",
            x_label="ECFP Tanimoto_avg (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )
        
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["2"])
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="ECFP Tanimoto_max (Coverage)",
            x_label="ECFP Tanimoto_max (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )
    
    st.write(f"radius={st.session_state["ecfp_settings"][0]}")
    st.write(f"fpSize={st.session_state["ecfp_settings"][1]}")

    st.divider()


    st.header(explanation[st.session_state["language"]]["common"]["6"])
    st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["3"])

    with st.spinner(explanation[st.session_state["language"]]["common"]["5"]):
        for i in range(st.session_state["n_conditions"]):
            st.subheader(st.session_state["label"][i])
            with st.container(border=True):

                st.write(
                    f"{explanation[st.session_state['language']]['ecfp_tanimoto']['5']}: {li_tanimoto_max[i]}"
                    )
            
                for index in li_max_indices[i]:
                    inchi1 = df1.iat[index[0], 0]
                    inchi2 = df1.iat[index[1], 0]

                    st.divider()
                    col1, col2 = st.columns(2)

                    with col1:
                        st.image(Draw.MolToImage(
                            Chem.MolFromInchi(inchi1), size=(600,300)
                            ))
                    with col2:
                        st.image(Draw.MolToImage(
                            Chem.MolFromInchi(inchi2), size=(600,300)
                            ))

                    st.caption(inchi1)
                    st.caption(inchi2)


    st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["4"])

    with st.spinner(explanation[st.session_state["language"]]["common"]["5"]):
        for i in range(st.session_state["n_conditions"]):
            st.subheader(st.session_state["label"][i])
            with st.container(border=True):

                st.write(
                    f"{explanation[st.session_state['language']]['ecfp_tanimoto']['5']}: {li_tanimoto_min[i]}"
                    )
            
                for index in li_min_indices[i]:
                    inchi1 = df1.iat[index[0], 0]
                    inchi2 = df1.iat[index[1], 0]

                    st.divider()
                    col1, col2 = st.columns(2)

                    with col1:
                        st.image(Draw.MolToImage(
                            Chem.MolFromInchi(inchi1), size=(600,300)
                            ))
                    with col2:
                        st.image(Draw.MolToImage(
                            Chem.MolFromInchi(inchi2), size=(600,300)
                            ))

                    st.caption(inchi1)
                    st.caption(inchi2)


    st.divider()
    
    st.subheader("ECFP Tanimoto")
    st.dataframe(df_result, hide_index=True)
