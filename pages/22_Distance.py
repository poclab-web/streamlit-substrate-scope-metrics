import io
import json

from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
import numpy as np
import pandas as pd
import streamlit as st
from streamlit_bokeh import streamlit_bokeh

from ssm import make_svg, metrics_distance
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="Distance")
st.title("Distance")

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

        li_distance_avg = []
        li_distance_max = []
        li_details = []
        for condition in st.session_state["conditions"]:
            df2a = df1[df1[condition] == 1].iloc[:, 0]
            df2b = df1[df1[condition] == 1].iloc[:, (st.session_state["n_conditions"] + 1):]
            df2 = pd.concat([df2a, df2b], axis=1)
            df2.reset_index(inplace=True, drop=True)

            ss = metrics_distance(df2)
            li_distance_avg.append(ss[2])
            li_distance_max.append(ss[4])
            li_details.append(ss[6])
            
        df2a = df1.iloc[:, 0]
        df2b = df1.iloc[:, (st.session_state["n_conditions"] + 1):]
        df2 = pd.concat([df2a, df2b], axis=1)

        ss = metrics_distance(df2)
        li_distance_avg.append(ss[2])
        li_distance_max.append(ss[4])


    df_result = pd.DataFrame(
        {"Condition": st.session_state["label"] + [f"Dataset ({len(df1)})"],
         "Distance_avg": li_distance_avg,
         "Distance_max": li_distance_max}
        )
    
    df_coverage = pd.DataFrame(
        {"Condition": st.session_state["label"],
         "Distance_avg (Coverage)": np.array(li_distance_avg[:-1]) / li_distance_avg[-1] * 100,
         "Distance_max (Coverage)": np.array(li_distance_max[:-1]) / li_distance_max[-1] * 100}
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
            y="Distance_avg",
            y_label="",
            horizontal=True,
            sort=False
            )
        
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["2"])
        st.bar_chart(
            df_result.iloc[:-1],
            x="Condition",
            y="Distance_max",
            y_label="",
            horizontal=True,
            sort=False
            )
        
    with tab2:
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["1"])
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="Distance_avg (Coverage)",
            x_label="Distance_avg (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )
        
        st.subheader(explanation[st.session_state["language"]]["ecfp_tanimoto"]["2"])
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="Distance_max (Coverage)",
            x_label="Distance_max (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )

    st.divider()


    st.header(explanation[st.session_state["language"]]["common"]["6"])
    st.subheader(explanation[st.session_state["language"]]["distance"]["1"])

    with st.spinner(explanation[st.session_state["language"]]["common"]["5"]):

        p = figure(
            x_range=(0, (df_result["Distance_max"].nlargest(2).min() + 1)),
            y_range=st.session_state["label"][::-1],
            height=300,
            x_axis_label="Distance",
            tools="hover, pan, box_zoom, reset, wheel_zoom"
            )

        for i in range(len(li_details)):
            df3 = li_details[i]
            df3 = df3[[df3.columns[0], (len(df3) - 1)]].iloc[:-1]
            df3.rename(columns={df3.columns[1]: "distance"}, inplace=True)

            df3["svg"] = df3.iloc[:, 0].map(make_svg)
            df3["condition"] = st.session_state["label"][i]

            p.scatter(
                x="distance",
                y="condition",
                source=ColumnDataSource(df3),
                size=7,
                line_color=None,
                fill_color="#21E6C1",
                hover_color="#FFF100"
            )

        p.xgrid.grid_line_color = None
        p.select_one(HoverTool).tooltips = """
        <div>
          @svg{safe}
        </div>
        <div>
          <span style="font-size: 20px; color: #1F4287;">Distance:</span>
          <span style="font-size: 20px; ">@distance</span>
        </div>
        """

    streamlit_bokeh(p)


    st.divider()
    
    st.subheader("Distance")
    st.dataframe(df_result, hide_index=True)
