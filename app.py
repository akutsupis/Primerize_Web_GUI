import streamlit as st
import os

# Import the tab UI modules
from tabs import tab1_baseline
from tabs import tab2_mapping
from tabs import tab3_structure
from tabs import tab4_custom
from tabs import tab5_protocols

# -----------------------------------------------------------------------------
# APPLICATION CONFIGURATION & SECURITY GUARDRAILS
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Primerize",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Reduce default blank space at the top of the page and enforce permanent scrollbar to prevent UI shifting
st.markdown(
    """
    <style>
    .block-container {
        padding-top: 3.5rem;
    }
    /* Enforce a permanent vertical scrollbar to prevent layout shifts on short tabs */
    html {
        overflow-y: scroll !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# -----------------------------------------------------------------------------
# HEADER & APPLICATION BRANDING
# -----------------------------------------------------------------------------
st.markdown(
    """
    **Primerize (previously named NA_thermo)** is a Python package for primer design and nucleic acid thermodynamics, developed by the Das Lab at Stanford University for high-throughput RNA synthesis and design.
    
    This website lets you run the Primerize algorithms through an intuitive graphical interface, with no coding required. You can design primers for simple 1D assembly, 2D chemical mapping libraries, and 3D structure-guided mutations, all in one place.
    
    *The original Primerize web server was decommissioned in May 2026. This website provides a graphical interface to the native Python backend.*
    
    **[Official Primerize Documentation & Tutorials](https://ribokit.github.io/Primerize/)**
    """
)
st.write("---")

# -----------------------------------------------------------------------------
# SIDEBAR PARAMETER SELECTION ENGINE (Shared Baseline Rules)
# -----------------------------------------------------------------------------
col1, col2, col3 = st.sidebar.columns([1, 2, 1])
with col2:
    st.image("logo_primerize.png", width='stretch')

st.sidebar.markdown("<br>", unsafe_allow_html=True)
st.sidebar.header("Baseline 1D Parameters")
sb_prefix = st.sidebar.text_input("Construct Prefix/Name", value="my_rna_construct")
sb_t7_check = st.sidebar.checkbox("Check for T7 promoter sequence", value=True)
sb_min_tm = st.sidebar.slider("Minimum Overlap Tm (°C)", min_value=45.0, max_value=85.0, value=60.0, step=0.5)

sb_limit_primers = st.sidebar.checkbox("Enforce Exact Primer Count Limit", value=False)
sb_num_primers = None 
if sb_limit_primers:
    sb_num_primers = st.sidebar.number_input("Exact Number of Primers (Even integers only)", min_value=2, max_value=100, value=8, step=2)

sb_min_len = st.sidebar.number_input("Minimum Primer Length (nt)", min_value=10, max_value=50, value=15)
sb_max_len = st.sidebar.number_input("Maximum Primer Length (nt)", min_value=40, max_value=120, value=60)

st.sidebar.info("To maintain stability and performance, sequences are limited to a maximum length of 350 nt. The Das Lab has tested the Primerize algorithm with sequences up to 300 nt.")

# -----------------------------------------------------------------------------
# MULTI-MODE TABS FOR ALL REPOSITORY FEATURES
# -----------------------------------------------------------------------------
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "1D Simple Assembly", 
    "2D Chemical Mapping", 
    "3D Structure Mutations", 
    "Custom Construct Mutants",
    "Protocols & Ordering"
])

# Route execution to modularized components
with tab1:
    tab1_baseline.render(sb_prefix, sb_t7_check, sb_min_tm, sb_num_primers, sb_min_len, sb_max_len)

with tab2:
    tab2_mapping.render()
    
with tab3:
    tab3_structure.render()
    
with tab4:
    tab4_custom.render()
    
with tab5:
    tab5_protocols.render()

# -----------------------------------------------------------------------------
# CITATION PANEL
# -----------------------------------------------------------------------------
st.write("---")
with st.expander("Academic Citation Acknowledgments"):
    st.markdown(
        """
        **This project is distributed under the MIT license.**

        Please cite the use of Primerize with the following foundational publications:
        - **Primerize-2D:** Tian, S., and Das, R. (2017) [Automated primer design for RNA multidimensional chemical mapping](https://academic.oup.com/bioinformatics/article-abstract/33/9/1405/2801460/Primerize-2D-automated-primer-design-for-RNA). *Bioinformatics* 33(9): 1405-1406.
        - **Primerize:** Tian, S., Yesselman, J.D., Cordero, P., and Das, R. (2015) [Primerize: automated primer assembly for transcribing interesting RNAs](http://nar.oxfordjournals.org/content/43/W1/W522). *Nucleic Acids Research* 43(W1): W522-W526.

        **Related Publications on the Mutate-Map-Rescue Pipeline:**
        - Tian, S., and Das, R. (2016) RNA structure through multidimensional chemical mapping. *Quarterly Review of Biophysics* 49(e7): 1-30.
        - Tian, S., Cordero, P., Kladwang, W., and Das, R. (2014) High-throughput mutate-map-rescue evaluates SHAPE-directed RNA structure and uncovers excited states. *RNA* 20(11): 1815-1826.
        - Kladwang, W., VanLang, C.C., Cordero P., and Das, R. (2011) A two-dimensional mutate-and-map strategy for non-coding RNA structure. *Nature Chemistry* 3: 954-962.
        """
    )