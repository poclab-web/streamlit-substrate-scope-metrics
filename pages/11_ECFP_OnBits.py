import io
import json

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, rdFingerprintGenerator
import streamlit as st

from ssm import generate_ecfp, on_bits
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="ECFP OnBits")
st.title("ECFP OnBits")

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
        df_ecfp = pd.DataFrame(np.array(st.session_state["ecfp"]))

        li_on_bits = []
        li_bit_list = []
        li_df2_ecfp = []
        for condition in st.session_state["conditions"]:
            df2_ecfp = df_ecfp[df1[condition] == 1].copy()

            ss = on_bits(df2_ecfp)
            li_on_bits.append(ss[0])
            li_bit_list.append(ss[1])
            li_df2_ecfp.append(ss[2])

        ss = on_bits(df_ecfp)
        li_on_bits.append(ss[0])


    df_result = pd.DataFrame(
        {"Condition": st.session_state["label"] + [f"Dataset ({len(df1)})"],
         "ECFP OnBits": li_on_bits}
        )
    
    df_coverage = pd.DataFrame(
        {"Condition": st.session_state["label"],
         "ECFP OnBits (Coverage)": np.array(li_on_bits[:-1]) / li_on_bits[-1] * 100}
        )
    
    
    tab1, tab2 = st.tabs([
        explanation[st.session_state["language"]]["common"]["7"],
        explanation[st.session_state["language"]]["common"]["8"]
    ])

    with tab1:
        st.bar_chart(
            df_result.iloc[:-1],
            x="Condition",
            y="ECFP OnBits",
            y_label="",
            horizontal=True,
            sort=False
            )
        
    with tab2:
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="ECFP OnBits (Coverage)",
            x_label="ECFP OnBits (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )
    
    st.write(f"radius={st.session_state["ecfp_settings"][0]}")
    st.write(f"fpSize={st.session_state["ecfp_settings"][1]}")

    st.divider()


    st.header(explanation[st.session_state["language"]]["common"]["6"])
    st.subheader(explanation[st.session_state["language"]]["ecfp_onbits"]["1"])

    with st.spinner(explanation[st.session_state["language"]]["common"]["5"]):
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(
            radius=st.session_state["ecfp_settings"][0],
            fpSize=st.session_state["ecfp_settings"][1]
            )
        additional_output = rdFingerprintGenerator.AdditionalOutput()

        for i in range(st.session_state["n_conditions"]):
            st.subheader(st.session_state["label"][i])
            
            draw_list = []
            for j in li_bit_list[i]:
                for k in li_df2_ecfp[i].index[:-1]:

                    if li_df2_ecfp[i].at[k, j] == 1:
                        mol = Chem.MolFromInchi(df1.iat[k, 0])
                        additional_output.CollectBitInfoMap()
                        mfpgen.GetFingerprint(mol, additionalOutput=additional_output)

                        draw_list.append((mol, j, additional_output.GetBitInfoMap()))
                        break
            
            img = Draw.DrawMorganBits(
                draw_list,
                molsPerRow=10,
                legends=[str(x) for x in li_bit_list[i]],
                useSVG=True
                )
            
            st.html(f"<style>.st-key-white{i}""{background-color: #FFFFFF}</style>")
            
            with st.container(key=f"white{i}"):
                st.image(img, width="stretch")


    st.divider()

    st.subheader("ECFP OnBits")
    st.dataframe(df_result, hide_index=True)
