import streamlit as st
import io
import sys
import os
import glob
import re
from utils import cached_design_custom, convert_ansi_to_html, strip_ansi, generate_advanced_files, generate_plate_svgs

def render():
    st.header("Custom Targeted Expression Library Controls")
    if 'active_job_1d' not in st.session_state:
        st.warning("Prerequisite Required: Establish a primary baseline assembly model in Tab 1 to run customized mutation operations.")
    else:
        t4_offset = st.number_input("Sequence Numbering Offset", value=0, key="t4_off")
        custom_mut_input = st.text_input("Explicit Targeted Mutants List", value="T120C, G119A;T120C", help="Separate variants using commas, semicolons, or whitespace blocks.")
        
        if st.button("Compile Tailored Custom Plate Layout", type="primary"):
            st.session_state['run_4d'] = True
            st.session_state['t4_inputs'] = (t4_offset, custom_mut_input)
            
        if st.session_state.get('run_4d', False):
            t4_offset, custom_mut_input = st.session_state['t4_inputs']
            with st.spinner("Assembling structural mutant mapping arrays..."):
                try:
                    job_custom, error_log, plates_data, zip_bytes = cached_design_custom(st.session_state['active_job_1d'], t4_offset, custom_mut_input)
                    if job_custom.is_success:
                        st.success("Custom Target Variant Mapping Completed.")
                        buf = io.StringIO()
                        sys.stdout = buf
                        print(job_custom)
                        sys.stdout = sys.__stdout__
                        output_text_custom = buf.getvalue()
                        
                        clean_html_output_custom = convert_ansi_to_html(output_text_custom)
                        st.markdown(
                            f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 18px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d; letter-spacing: 0.5px;">{clean_html_output_custom}</div>', 
                            unsafe_allow_html=True
                        )
                        st.markdown("<br>", unsafe_allow_html=True)
                        
                        col_dl1, col_dl2 = st.columns(2)
                        with col_dl1:
                            st.download_button("Download Custom Plate Specification File (.txt)", data=strip_ansi(output_text_custom), file_name=f"{st.session_state['active_job_1d'].name}_custom_plate.txt")
                            
                        with col_dl2:
                            xls_name = f"{st.session_state['active_job_1d'].name}_custom"
                            job_custom.save('table', path='.', name=xls_name)
                            xls_files = glob.glob(f"{xls_name}_plate_*.xls")
                            for i, xls_path in enumerate(xls_files):
                                with open(xls_path, "rb") as f:
                                    xls_data = f.read()
                                st.download_button(f"Download Plate {i+1} (.xls)", data=xls_data, file_name=xls_path, mime="application/vnd.ms-excel", key=f"dl_4d_{i}")
                                os.remove(xls_path)
                                
                        with st.expander("Supplementary Data Files"):
                            st.caption("Supplementary output files for downstream analysis and sequence alignment pipelines.")
                            c1, c2, c3 = st.columns(3)
                            constructs_txt, assembly_txt, _ = generate_advanced_files(job_custom)
                            with c1:
                                st.download_button("Download Constructs Map (.txt)", data=constructs_txt, file_name=f"{job_custom.name}_constructs.txt", key="t4_dl_const")
                            with c2:
                                st.download_button("Download Assembly DNA (.txt)", data=assembly_txt, file_name=f"{job_custom.name}_assembly.txt", key="t4_dl_assem")

                        if plates_data:
                            st.write("---")
                            st.subheader("Plate Layout")
                            cols = st.columns(4)
                            for idx, p_data in enumerate(plates_data):
                                with cols[idx % 4]:
                                    st.markdown(f"<p style='text-align: center; font-weight: bold;'>{p_data['plate_name']}</p>", unsafe_allow_html=True)
                                    st.markdown(p_data['svg_string'], unsafe_allow_html=True)
                            st.markdown("<br>", unsafe_allow_html=True)
                            st.download_button("Download Plate Graphics (ZIP)", data=zip_bytes, file_name=f"{job_custom.name}_plates.zip", mime="application/zip", key="t4_dl_plates_zip")
                                
                    else:
                        clean_err = strip_ansi(error_log).strip()
                        val_err_match = re.search(r'ValueError:\s*(?:ERROR:\s*)?(.*)', clean_err)
                        if val_err_match:
                            parsed_err = val_err_match.group(1).strip()
                            st.error(f"Error: Custom Construction Intercept. The target variations could not map to the active sequence range safely.\n\n**{parsed_err}**")
                        else:
                            st.error(f"Error: Custom Construction Intercept. The target variations could not map to the active sequence range safely.\n\n**Backend Output:**\n`{clean_err}`")
                except Exception as e:
                    st.error(f"Error Compiling Custom Plate Array: {str(e)}")
