import streamlit as st
import io
import sys
import os
import glob
from utils import cached_design_2d, convert_ansi_to_html, strip_ansi, generate_advanced_files

def render():
    st.header("Automated High-Throughput Mutagenesis (2D Libraries)")
    if 'active_job_1d' not in st.session_state:
        st.warning("Prerequisite Required: You must generate a valid 1D baseline design in Tab 1 before building a secondary 2D library layout.")
    else:
        baseline = st.session_state['active_job_1d']
        st.info(f"Loaded Active Baseline Sequence: **{baseline.name}** ({len(baseline.sequence)} nt)")
        
        seq_len = len(baseline.sequence)
        
        col1, col2, col3 = st.columns(3)
        with col1:
            t2_offset = st.number_input("Sequence Numbering Offset", value=0, help="Sequence numbering offset, which is one minus the final number of the first nucleotide.")
        with col2:
            t2_min_mut = st.number_input("Mutation Start Index", min_value=1, max_value=seq_len, value=1)
        with col3:
            t2_max_mut = st.number_input("Mutation End Index", min_value=1, max_value=seq_len, value=seq_len)
            
        t2_lib = st.selectbox("Target Library Assembly Strategy (which_lib)", options=[1, 2, 3], index=0)
        
        if st.button("Generate 2D Mapping Library", type="primary"):
            st.session_state['run_2d'] = True
            st.session_state['t2_inputs'] = (t2_offset, t2_min_mut, t2_max_mut, t2_lib)
            
        if st.session_state.get('run_2d', False):
            t2_offset, t2_min_mut, t2_max_mut, t2_lib = st.session_state['t2_inputs']
            with st.spinner("Processing mutation sequence matrices..."):
                try:
                    job_2d = cached_design_2d(baseline, t2_offset, t2_min_mut, t2_max_mut, t2_lib)
                    if job_2d.is_success:
                        st.success("2D Plate Mapping Arrays Generated Successfully.")
                        
                        buf = io.StringIO()
                        sys.stdout = buf
                        print(job_2d)
                        sys.stdout = sys.__stdout__
                        output_text_2d = buf.getvalue()

                        clean_html_output_2d = convert_ansi_to_html(output_text_2d)

                        st.markdown(
                            f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 18px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d; letter-spacing: 0.5px;">{clean_html_output_2d}</div>', 
                            unsafe_allow_html=True
                        )
                        st.markdown("<br>", unsafe_allow_html=True)
                        
                        col_dl1, col_dl2 = st.columns(2)
                        with col_dl1:
                            st.download_button("Download 2D Plate Specification File (.txt)", data=strip_ansi(output_text_2d), file_name=f"{baseline.name}_2D_library.txt")
                            
                        with col_dl2:
                            xls_name = f"{baseline.name}_2D_library"
                            job_2d.save('table', path='.', name=xls_name)
                            xls_files = glob.glob(f"{xls_name}_plate_*.xls")
                            for i, xls_path in enumerate(xls_files):
                                with open(xls_path, "rb") as f:
                                    xls_data = f.read()
                                st.download_button(f"Download Plate {i+1} (.xls)", data=xls_data, file_name=xls_path, mime="application/vnd.ms-excel", key=f"dl_2d_{i}")
                                os.remove(xls_path)
                                
                        with st.expander("Supplementary Data Files"):
                            st.caption("Supplementary output files for downstream analysis and sequence alignment pipelines.")
                            c1, c2, c3 = st.columns(3)
                            constructs_txt, assembly_txt, _ = generate_advanced_files(job_2d)
                            with c1:
                                st.download_button("Download Constructs Map (.txt)", data=constructs_txt, file_name=f"{job_2d.name}_constructs.txt", key="t2_dl_const")
                            with c2:
                                st.download_button("Download Assembly DNA (.txt)", data=assembly_txt, file_name=f"{job_2d.name}_assembly.txt", key="t2_dl_assem")
                                
                    else:
                        st.error("Error: 2D Processing Failure. No valid plate layout found for this specific mutation range.")
                except Exception as e:
                    st.error(f"Error: {str(e)}")
