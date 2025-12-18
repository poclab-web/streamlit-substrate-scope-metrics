import io
import json

from bokeh.models import ColumnDataSource, HoverTool, Legend, LegendItem
from bokeh.plotting import figure
import numpy as np
import pandas as pd
import streamlit as st
from streamlit_bokeh import streamlit_bokeh

from ssm import inchi_to_svg, metrics_ConvexHull, pca_2d
from utils.sidebar import sidebar


sidebar()
st.set_page_config(page_title="ConvexHull")
st.title("ConvexHull")

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

        if st.session_state["pca"] is None:
            df_pca1, ratio = pca_2d(
                df1.iloc[:, (st.session_state["n_conditions"] + 1):]
                )
            df_pca2 = pd.concat(
                [df1.iloc[:, :(st.session_state["n_conditions"] + 1)], df_pca1],
                axis=1
                )
            
            st.session_state["pca"] = df_pca2.to_csv(index=False)
            st.session_state["pca_ratio"] = ratio


        df_pca = pd.read_csv(io.StringIO(st.session_state["pca"]))
        df_pca["svg"] = inchi_to_svg(df_pca.iloc[:, 0])
        
        li_convexhull = []
        li_bokeh = []
        for condition in st.session_state["conditions"]:
            df2_pca = df_pca[df_pca[condition] == 1].copy()
            array2 = df2_pca.loc[:, "PC1":"PC2"].to_numpy()

            area, vertices = metrics_ConvexHull(array2)
            li_convexhull.append(area)


            df_u = df_pca[df_pca[condition] == 0]
            df_n = df_pca[df_pca[condition] == -1]

            p = figure(
                x_axis_label="PC1",
                y_axis_label="PC2",
                tools="hover, pan, box_zoom, reset, wheel_zoom"
                )
            item_list = []

            r4 = p.patch(
                x=array2[vertices].T[0],
                y=array2[vertices].T[1],
                color="#21E6C1",
                alpha=0.3
                )
            item_list.insert(0, LegendItem(label="ConvexHull", renderers=[r4]))
            
            if len(df_u) > 0:
                r3 = p.scatter(
                    x="PC1",
                    y="PC2",
                    source=ColumnDataSource(df_u[["PC1", "PC2", "svg"]]),
                    size=3,
                    line_color=None,
                    fill_color="#C8C8CB",
                    hover_color="#FFF100",
                    alpha=0.5
                    )
                item_list.insert(0, LegendItem(label="Unlabeled", renderers=[r3]))
            
            if len(df_n) > 0:
                r2 = p.scatter(
                    x="PC1",
                    y="PC2",
                    source=ColumnDataSource(df_n[["PC1", "PC2", "svg"]]),
                    size=5,
                    line_color=None,
                    fill_color="#FFFFFF",
                    hover_color="#FFF100"
                    )
                item_list.insert(0, LegendItem(label="Negative", renderers=[r2]))
            
            r1 = p.scatter(
                x="PC1",
                y="PC2",
                source=ColumnDataSource(df2_pca[["PC1", "PC2", "svg"]]),
                size=7,
                line_color=None,
                fill_color="#21E6C1",
                hover_color="#FFF100"
                )
            item_list.insert(0, LegendItem(label="Positive", renderers=[r1]))
            
            p.add_layout(Legend(items=item_list))
            p.select_one(HoverTool).tooltips = "<div>@svg{safe}</div>"
            
            li_bokeh.append(p)


        array = df_pca.loc[:, "PC1":"PC2"].to_numpy()
        area, _ = metrics_ConvexHull(array)
        li_convexhull.append(area)
        
        
    df_result = pd.DataFrame(
        {"Condition": st.session_state["label"] + [f"Dataset ({len(df_pca)})"],
        "ConvexHull": li_convexhull}
        )
    
    df_coverage = pd.DataFrame(
        {"Condition": st.session_state["label"],
        "ConvexHull (Coverage)": np.array(li_convexhull[:-1]) / li_convexhull[-1] * 100}
        )


    tab1, tab2 = st.tabs([
        explanation[st.session_state["language"]]["common"]["7"],
        explanation[st.session_state["language"]]["common"]["8"]
    ])

    with tab1:
        st.bar_chart(
            df_result.iloc[:-1],
            x="Condition",
            y="ConvexHull",
            y_label="",
            horizontal=True,
            sort=False
            )
        
    with tab2:
        st.bar_chart(
            df_coverage,
            x="Condition",
            y="ConvexHull (Coverage)",
            x_label="ConvexHull (Coverage) [%]",
            y_label="",
            horizontal=True,
            sort=False
            )

    st.divider()


    st.header(explanation[st.session_state["language"]]["common"]["6"])

    with st.spinner(explanation[st.session_state["language"]]["common"]["5"]):
        st.subheader(explanation[st.session_state["language"]]["convexhull"]["1"])
        st.write(f"PC1: {st.session_state['pca_ratio'][0] * 100} %")
        st.write(f"PC2: {st.session_state['pca_ratio'][1] * 100} %")

        for i in range(st.session_state["n_conditions"]):
            st.subheader(st.session_state["label"][i])
            streamlit_bokeh(li_bokeh[i])

        st.subheader(explanation[st.session_state["language"]]["convexhull"]["3"])
        df_corr = pd.concat(
            [df_pca.loc[:, "PC1":"PC2"],
             df1.iloc[:, (st.session_state["n_conditions"] + 1):]],
            axis=1
        ).corr()
        st.dataframe(df_corr.iloc[2:, 0:2])

    st.subheader(explanation[st.session_state["language"]]["convexhull"]["2"])
    df_pca.rename(columns={df_pca.columns[0]: "InChI"}, inplace=True)
    st.dataframe(df_pca.iloc[:, :-1], hide_index=True)
    
    
    st.divider()
    
    st.subheader("ConvexHull")
    st.dataframe(df_result, hide_index=True)
