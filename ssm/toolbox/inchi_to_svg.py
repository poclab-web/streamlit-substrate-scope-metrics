from multiprocessing import Pool, cpu_count

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import streamlit as st


def make_svg(inchi):
    mol = Chem.MolFromInchi(inchi)
    svg = rdMolDraw2D.MolToSVG(mol, width=300, height=150)

    return svg


@st.cache_data
def inchi_to_svg(inchi_list):
    with Pool(cpu_count()) as p:
        svg_list = p.map(make_svg, inchi_list)

    return svg_list
