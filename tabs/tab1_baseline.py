import streamlit as st
import pandas as pd
import io
import sys
from utils import cached_design_1d, convert_ansi_to_html, strip_ansi

def render(sb_prefix, sb_min_tm, sb_num_primers, sb_min_len, sb_max_len):
    st.header("Simple PCR Assembly Setup")
    raw_seq = st.text_area("Target Sequence (DNA or RNA)", placeholder="PASTE YOUR SEQUENCE HERE...", height=180, key="t1_seq")
    
    if st.button("Generate Baseline 1D Assembly", type="primary", key="btn_1d"):
        st.session_state['run_1d'] = True
        st.session_state['raw_seq_1d'] = raw_seq
        
    if st.session_state.get('run_1d', False):
        cleaned_seq = "".join(st.session_state['raw_seq_1d'].split()).upper().replace('U', 'T')
        
        if not cleaned_seq:
            st.error("Submission Failure: The sequence text boundary cannot be blank.")
        elif len(cleaned_seq) > 350:
            st.error(f"Length Alert: The sequence is {len(cleaned_seq)} nt long. To keep free tier execution fast and stable, sequences are limited to 350 nt.")
        else:
            with st.spinner("Calculating optimal primer boundaries..."):
                try:
                    job_1d = cached_design_1d(cleaned_seq, sb_min_tm, sb_num_primers, sb_min_len, sb_max_len, sb_prefix)
                    
                    if job_1d.is_success:
                        st.success("Assembly primer optimization completed successfully.")
                        st.session_state['active_job_1d'] = job_1d
                        
                        st.subheader("Primary Assembly Primer Set")
                        primer_records = [{"Primer Index": i + 1, "Oligo Sequence (5' → 3')": s} for i, s in enumerate(job_1d.primer_set)]
                        st.dataframe(pd.DataFrame(primer_records), width="stretch")
                        
                        st.subheader("Full Backend Log")
                        st.caption("Raw output engine execution log.")
                        
                        buffer = io.StringIO()
                        sys.stdout = buffer
                        print(job_1d)
                        sys.stdout = sys.__stdout__
                        output_text = buffer.getvalue()
                        
                        clean_html_output = convert_ansi_to_html(output_text)

                        st.markdown(
                            f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 18px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d; letter-spacing: 0.5px;">{clean_html_output}</div>', 
                            unsafe_allow_html=True
                        )
                        st.markdown("<br>", unsafe_allow_html=True)
                            
                        st.download_button("Download Design Matrix Text File", data=strip_ansi(output_text), file_name=f"{job_1d.name}_1D_design.txt")
                        
                        st.write("---")
                        st.subheader("IDT Bulk Ordering Block")
                        st.info("Copy the text block below and paste it directly into the IDT Bulk Input portal. Select 'Lab Ready' for normalization.")
                        
                        idt_lines = [
                            '#',
                            '',
                            '------/* IDT USER: for primer ordering, copy and paste to Bulk Input */------',
                            '------/* START */------'
                        ]
                        for i in range(len(job_1d.primer_set)):
                            suffix = 'FR'[i % 2]
                            idt_lines.append(f'{job_1d.name}-{i + 1}{suffix}\t{job_1d.primer_set[i]}\t\t25nm\tSTD')
                        idt_lines.extend([
                            '------/* END */------',
                            '------/* NOTE: use "Lab Ready" for "Normalization" */------'
                        ])
                        idt_text = '\n'.join(idt_lines)
                        st.code(idt_text, language='text')
                        
                    else:
                        st.error("Error: Primerize Engine Failure. No valid primer boundaries could be evaluated under these parameters.")
                except Exception as e:
                    st.error(f"Execution Error Intercepted: {str(e)}")
