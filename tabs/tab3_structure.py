import streamlit as st
import io
import sys
import os
import glob
import pandas as pd
import io
import sys
import primerize
from utils import cached_design_3d, convert_ansi_to_html, strip_ansi, generate_advanced_files

def render():
    st.header("3D Structure Mutations")
    st.caption("Optimize assembly primer sets for structural targets by integrating thermodynamic calculations with dot-bracket folding models.")

    if 'active_job_1d' not in st.session_state:
        st.warning(
            "**Prerequisite Required:** You must run a 1D Baseline calculation in Tab 1 "
            "to unlock the 3D structure-guided modules. This populates the foundational "
            "primer landscape required for mutation mapping."
        )
    else:
        job_1d = st.session_state['active_job_1d']
        st.info(f"Active Baseline: {job_1d.name} ({len(job_1d.sequence)} nt)")
        
        seq_3d = st.text_area(
            "Target Sequence (Synchronized from Tab 1)", 
            value=job_1d.sequence,
            disabled=True,
            help="This sequence is locked to match your active 1D baseline project.",
            height=120,
            key="t3_seq_input"
        )
        
        struct_3d = st.text_area(
            "Dot-Bracket Structure Map (e.g., (((...))) )", 
            placeholder="Paste matching secondary structure fold layout configuration here...",
            help="Enter a valid secondary structure string consisting of dots (.) and matching brackets ().",
            height=120,
            key="t3_struct_input"
        )
        
        if st.button("Generate 3D Mutation Design", type="primary", key="btn_3d"):
            st.session_state['run_3d'] = True
            st.session_state['t3_struct'] = struct_3d
            
        if st.session_state.get('run_3d', False):
            cleaned_struct = "".join(st.session_state['t3_struct'].split())
            
            if not cleaned_struct:
                st.error("Submission Failure: The Dot-Bracket Structure map cannot be blank.")
            elif len(job_1d.sequence) != len(cleaned_struct):
                st.error(
                    f"**Input Length Mismatch!**\n\n"
                    f"The 3D structural folding map must match your target sequence length character-for-character:\n"
                    f"* Your target sequence is **{len(job_1d.sequence)}** nucleotides long.\n"
                    f"* Your dot-bracket structure input is **{len(cleaned_struct)}** characters long.\n\n"
                    f"Please adjust your structure layout to align perfectly before running the pipeline."
                )
            else:
                with st.spinner("Calculating 3D structure-guided mutations and primer boundaries..."):
                    try:
                        job_3d, plates_data, zip_bytes = cached_design_3d(
                            _job_1d=job_1d, 
                            offset=0, 
                            structures=[cleaned_struct], 
                            n_mutations=2, 
                            which_lib=1, 
                            is_single=True, 
                            is_fillwt=True
                        )
                        
                        if job_3d.is_success:
                            st.success("3D Structural Synthesis Optimization Completed.")
                            
                            st.subheader("Structural Mutation Assembly Primers")
                            
                            hcol1, hcol2, hcol3 = st.columns([1, 1, 8])
                            with hcol1: st.markdown("**#**")
                            with hcol2: st.markdown("**Length**")
                            with hcol3: st.markdown("**Sequence** (To copy, hover over a row and click the icon on the far right)")
                            st.markdown("---")
                            
                            for i, primer in enumerate(job_3d.primer_set):
                                col1, col2, col3 = st.columns([1, 1, 8])
                                
                                with col1:
                                    badge = "F" if i % 2 == 0 else "R"
                                    color = "#00bcd4" if badge == "F" else "#ff5722"
                                    st.markdown(f"<span style='background-color:{color}; color:white; padding: 3px 8px; border-radius: 4px; font-weight: bold; font-family: monospace;'>{i+1} [{badge}]</span>", unsafe_allow_html=True)
                                with col2:
                                    st.markdown(f"**{len(primer)}**")
                                with col3:
                                    st.code(primer, language="text")
                            
                            buffer_3d = io.StringIO()
                            sys.stdout = buffer_3d
                            print(job_3d)
                            sys.stdout = sys.__stdout__
                            output_text_3d = buffer_3d.getvalue()
                            
                            clean_html_output_3d = convert_ansi_to_html(output_text_3d)
                            
                            st.subheader("3D Alignment Matrix Mapping")
                            st.markdown(
                                f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 15px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d;">{clean_html_output_3d}</div>', 
                                unsafe_allow_html=True
                            )
                            st.markdown("<br>", unsafe_allow_html=True)
                            
                            col_dl1, col_dl2 = st.columns(2)
                            with col_dl1:
                                st.download_button(
                                    "Download 3D Design Matrix (.txt)", 
                                    data=strip_ansi(output_text_3d), 
                                    file_name=f"{job_1d.name}_3D_structural_design.txt"
                                )
                                
                            with col_dl2:
                                xls_name = f"{job_1d.name}_3D_structural_design"
                                job_3d.save('table', path='.', name=xls_name)
                                xls_files = glob.glob(f"{xls_name}_plate_*.xls")
                                for i, xls_path in enumerate(xls_files):
                                    with open(xls_path, "rb") as f:
                                        xls_data = f.read()
                                    st.download_button(f"Download Plate {i+1} (.xls)", data=xls_data, file_name=xls_path, mime="application/vnd.ms-excel", key=f"dl_3d_{i}")
                                    os.remove(xls_path)
                                    
                            with st.expander("Supplementary Data Files"):
                                st.caption("Supplementary output files for downstream analysis and sequence alignment pipelines.")
                                c1, c2, c3 = st.columns(3)
                                constructs_txt, assembly_txt, structs_txt = generate_advanced_files(job_3d, include_structures=True, structures_list=[cleaned_struct])
                                with c1:
                                    st.download_button("Download Constructs Map (.txt)", data=constructs_txt, file_name=f"{job_3d.name}_constructs.txt", key="t3_dl_const")
                                with c2:
                                    st.download_button("Download Assembly DNA (.txt)", data=assembly_txt, file_name=f"{job_3d.name}_assembly.txt", key="t3_dl_assem")
                                with c3:
                                    st.download_button("Download Structures Map (.txt)", data=structs_txt, file_name=f"{job_3d.name}_structures.txt", key="t3_dl_struct")
                                    
                            # 8. Visual Plate Layouts
                            if plates_data:
                                st.write("---")
                                st.subheader("Plate Layout")
                                
                                cols = st.columns(4)
                                for idx, p_data in enumerate(plates_data):
                                    with cols[idx % 4]:
                                        st.markdown(f"<p style='text-align: center; font-weight: bold;'>{p_data['primer_name']}</p>", unsafe_allow_html=True)
                                        st.markdown(p_data['svg_string'], unsafe_allow_html=True)
                                        
                                st.markdown("<br>", unsafe_allow_html=True)
                                st.download_button("Download Plate Graphics (ZIP)", data=zip_bytes, file_name=f"{job_3d.name}_plates.zip", mime="application/zip", key="t3_dl_plates_zip")
                                    
                        else:
                            st.error("Error: Primerize Engine Failure. No valid primer combinations could satisfy the 3D folding constraints.")
                            
                    except Exception as e:
                        st.error(f"Execution Error Intercepted: {str(e)}")
