import streamlit as st
import pandas as pd
import io
import sys
from utils import cached_design_1d, convert_ansi_to_html, strip_ansi

def render(sb_prefix, sb_t7_check, sb_min_tm, sb_num_primers, sb_min_len, sb_max_len):
    st.header("Simple PCR Assembly Setup")
    
    st.markdown("""
    **Target Sequence (DNA or RNA):**
    - Please enter your sequence below: nucleotides only, no headers or comments.
    - Valid nucleotides are A, C, G, T, and U; and at least 60 nt long.
    - Flanking sequences (e.g. T7 promoter, buffering region, tail) should be included.
    """)
    raw_seq = st.text_area("Target Sequence (DNA or RNA)", placeholder="PASTE YOUR SEQUENCE HERE...", height=180, key="t1_seq", label_visibility="collapsed")
    
    if st.button("Primerize!", type="primary", key="btn_1d"):
        st.session_state['run_1d'] = True
        st.session_state['raw_seq_1d'] = raw_seq
        
    if st.session_state.get('run_1d', False):
        cleaned_seq = "".join(st.session_state['raw_seq_1d'].split()).upper().replace('U', 'T')
        
        if not cleaned_seq:
            st.error("Submission Failure: The sequence text boundary cannot be blank.")
            return

        # 1. T7 Promoter Logic
        t7_action_msg = "T7 Promoter Sequence: Feature disabled."
        if sb_t7_check:
            t7_seq = "TTCTAATACGACTCACTATAGGG"
            if not cleaned_seq.startswith(t7_seq):
                cleaned_seq = t7_seq + cleaned_seq
                t7_action_msg = "T7 Promoter Sequence: Feature enabled (uncheck the option to disable). T7 promoter sequence was missing, sequence has been added."
            else:
                t7_action_msg = "T7 Promoter Sequence: Feature enabled (uncheck the option to disable). T7 promoter sequence was present, no action was taken."
        
        # 2. Display Length & Run Engine
        if len(cleaned_seq) > 350:
            st.error(f"Length Alert: The sequence is {len(cleaned_seq)} nt long. To keep free tier execution fast and stable, sequences are limited to 350 nt.")
        else:
            with st.spinner("Calculating optimal primer boundaries..."):
                try:
                    job_1d = cached_design_1d(cleaned_seq, sb_min_tm, sb_num_primers, sb_min_len, sb_max_len, sb_prefix)
                    
                    if job_1d.is_success:
                        st.success(f"Assembly primer optimization completed successfully for **{len(cleaned_seq)} nt** sequence.")
                        st.session_state['active_job_1d'] = job_1d
                        
                        # 3. T7 Status Indicator
                        st.info(t7_action_msg)

                        # 4. Capture and parse Backend Log for Warnings
                        buffer = io.StringIO()
                        sys.stdout = buffer
                        print(job_1d)
                        sys.stdout = sys.__stdout__
                        raw_output = buffer.getvalue()
                        
                        warnings = []
                        clean_log_lines = []
                        import re
                        for line in raw_output.splitlines():
                            clean_line = strip_ansi(line).strip()
                            if clean_line.startswith("WARNING:"):
                                # Make 'WARNING:' bold
                                clean_line = clean_line.replace("WARNING:", "**WARNING:**")
                                # Inject Streamlit native markdown colors for simplicity and maintainability
                                # Escaping the inner brackets so Streamlit's markdown parser doesn't break
                                clean_line = re.sub(r'(\d+)\s+F\b', r":blue[\1 \[F\]]", clean_line)
                                clean_line = re.sub(r'(\d+)\s+R\b', r":red[\1 \[R\]]", clean_line)
                                warnings.append(clean_line)
                            else:
                                clean_log_lines.append(line)
                                
                        if warnings:
                            # Append the generic web-server UI warning if any specific mispriming warnings fired
                            warnings.append("**WARNING:** One-pot PCR assembly may fail due to mispriming; consider first assembling fragments in a preliminary PCR round (subpool).")
                            st.warning("\n\n".join(warnings))
                        
                        # 5. Visual Primer Set Table with Badges & Copy-to-Clipboard
                        st.subheader("Designed Primers")
                        
                        # Create visual column headers
                        hcol1, hcol2, hcol3 = st.columns([1, 1, 8])
                        with hcol1: st.markdown("**#**")
                        with hcol2: st.markdown("**Length**")
                        with hcol3: st.markdown("**Sequence** (To copy, hover over a row and click the icon on the far right)")
                        st.markdown("---")
                        
                        for i, primer in enumerate(job_1d.primer_set):
                            col1, col2, col3 = st.columns([1, 1, 8])
                            
                            with col1:
                                badge = "F" if i % 2 == 0 else "R"
                                color = "#00bcd4" if badge == "F" else "#ff5722"
                                st.markdown(f"<span style='background-color:{color}; color:white; padding: 3px 8px; border-radius: 4px; font-weight: bold; font-family: monospace;'>{i+1} [{badge}]</span>", unsafe_allow_html=True)
                            with col2:
                                st.markdown(f"**{len(primer)}**")
                            with col3:
                                st.code(primer, language="text")
                                
                        # 6. Assembly Scheme (Clean Log)
                        st.subheader("Assembly Scheme")
                        clean_html_output = convert_ansi_to_html("\n".join(clean_log_lines))
                        st.markdown(
                            f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 18px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d; letter-spacing: 0.5px;">{clean_html_output}</div>', 
                            unsafe_allow_html=True
                        )
                        st.markdown("<br>", unsafe_allow_html=True)
                        
                        dcol1, dcol2 = st.columns(2)
                        with dcol1:
                            st.download_button("Download Design Matrix Text File", data=strip_ansi(raw_output), file_name=f"{job_1d.name}_1D_design.txt")
                        with dcol2:
                            csv_data = "Primer Name,Sequence\n"
                            for i, p in enumerate(job_1d.primer_set):
                                suffix = 'FR'[i % 2]
                                csv_data += f"{job_1d.name}-{i + 1}{suffix},{p}\n"
                            st.download_button("Download Primers CSV (Tubes)", data=csv_data.encode('utf-8'), file_name=f"{job_1d.name}_primers.csv", mime="text/csv")
                        
                        # 7. Cleaned IDT Bulk Ordering Block
                        st.write("---")
                        st.subheader("IDT Bulk Ordering Block")
                        
                        idt_lines = [
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
