import streamlit as st
import pandas as pd
import io
import re
import sys
import os
import threading
from contextlib import redirect_stdout
from ansi2html import Ansi2HTMLConverter

# Ensure the local vendored primerize package takes import precedence
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import primerize

# Global thread lock to prevent Singleton race conditions in the primerize backend
primerize_lock = threading.Lock()

def strip_ansi(text):
    """Removes ANSI escape sequences from text for clean raw file downloads."""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return re.sub(ansi_escape, '', text)

def generate_advanced_files(job, include_structures=False, structures_list=None):
    # 1. Constructs Map
    try:
        lib_num = job.get('which_lib')
    except AttributeError:
        lib_num = 1
    prefix = f"Lib{lib_num}-" if lib_num else 'Lib1-'
    
    constructs_clean = ""
    if 'constructs' in job._data:
        constructs_text = job._data['constructs'].echo(prefix)
        constructs_clean = strip_ansi(constructs_text)
        
    # 2. Assembly Map
    assembly_clean = ""
    if 'assembly' in job._data:
        assembly_lines = job._data['assembly'].echo()
        assembly_lines += '\nPRIMERS             LENGTH    \tSEQUENCE\n'
        for i, primer in enumerate(job.primer_set):
            suffix = ' R' if i % 2 else ' F'
            name_str = f"{job.name}-{i+1}{suffix}"
            assembly_lines += f"{name_str.ljust(39)}{str(len(primer)).ljust(10)}{primer}\n"
        assembly_clean = strip_ansi(assembly_lines)
        
    # 3. Structures Map (Optional)
    structures_clean = ""
    if include_structures and structures_list:
        lines = list(structures_list)
        if 'warnings' in job._data and job._data['warnings']:
            lines.extend(['', 'WARNINGS:', ''])
            try:
                offset = job.get('offset')
            except AttributeError:
                offset = 0
            seq = job.sequence
            for pair in job._data['warnings']:
                lines.append(f"Mismatch in base-pair between {seq[pair[0] - 1]}{pair[0] - offset} and {seq[pair[1] - 1]}{pair[1] - offset}.")
        structures_clean = '\n'.join(lines)
        
    return constructs_clean, assembly_clean, structures_clean

# Helper function to convert ansi to HTML for Streamlit rendering
from ansi2html import Ansi2HTMLConverter

def convert_ansi_to_html(text):
    """Converts terminal ANSI codes to clean, secure, and properly bound inline HTML."""
    if not text:
        return ""
    
    # 1. Clean out legacy text file carriage components
    cleaned_text = text.replace("\r\n", "\n").replace("\r", "\n")
    
    # 2. Instantiate converter without full page wrapper headers or external stylesheets
    converter = Ansi2HTMLConverter(inline=True, linkify=False)
    html_core = converter.convert(cleaned_text, full=False)
    
    # 3. Strip out any stray tag literals that might confuse the Streamlit Markdown engine
    html_core = html_core.replace("</div>", "").replace("<div>", "")
    
    return html_core.strip()

def strip_ansi(text):
    """Removes ANSI escape sequences from text for clean raw file downloads."""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)

# -----------------------------------------------------------------------------
# APPLICATION CONFIGURATION & SECURITY GUARDRAILS
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Primerize",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Reduce default blank space at the top of the page without clipping headers
st.markdown(
    """
    <style>
    .block-container {
        padding-top: 3.5rem;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# Cache execution results to prevent redundant calculations across UI re-renders
@st.cache_data(max_entries=128, ttl=3600)
def cached_design_1d(seq, min_tm, num_primers, min_len, max_len, prefix):
    prm_1d = primerize.Primerize_1D
    
    # 1. Reset the singleton back to factory defaults to clear out old runs
    prm_1d.reset()
    
    # 2. Pass the parameters directly into the design function call execution layer.
    # Note: We must bypass the module's native .set() method here because it contains
    # a known architectural flaw where it fails to properly initialize variables in 
    # the singleton state, leading to runtime ReferenceErrors.
    with primerize_lock:
        return prm_1d.design(
            sequence=seq,
            MIN_TM=float(min_tm),
            NUM_PRIMERS=int(num_primers) if num_primers else None,
            MIN_LENGTH=int(min_len),
            MAX_LENGTH=int(max_len),
            prefix=str(prefix)
        )

@st.cache_data(max_entries=64, ttl=1800)
def cached_design_2d(_job_1d, offset, min_mut, max_mut, which_lib):
    mut_range = list(range(int(min_mut), int(max_mut) + 1))
    with primerize_lock:
        return primerize.Primerize_2D.design(_job_1d, offset=int(offset), which_muts=mut_range, which_lib=int(which_lib))

@st.cache_data(max_entries=64, ttl=1800)
def cached_design_3d(_job_1d, offset, structures, n_mutations, which_lib, is_single, is_fillwt):
    with primerize_lock:
        return primerize.Primerize_3D.design(
            _job_1d, 
            offset=int(offset), 
            structures=structures, 
            N_mutations=int(n_mutations), 
            which_lib=int(which_lib), 
            is_single=is_single, 
            is_fillWT=is_fillwt
        )

@st.cache_data(max_entries=64, ttl=1800)
def cached_design_custom(_job_1d, offset, raw_mutation_string):
    mut_list = primerize.Construct_List()
    # Clean and split string input (supports comma, space, or semicolon separated entries)
    parsed_muts = [m.strip() for m in re.split(r'[,\s;]+', raw_mutation_string) if m.strip()]
    
    # Push individual targeted construct arrays to the execution engine
    mut_list.push(parsed_muts)
    
    f = io.StringIO()
    with primerize_lock:
        with redirect_stdout(f):
            job = primerize.Primerize_Custom.design(_job_1d, offset=int(offset), mut_list=mut_list)
            
    return job, f.getvalue()

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
# Use columns to shrink and vertically center the logo in the sidebar
col1, col2, col3 = st.sidebar.columns([1, 2, 1])
with col2:
    st.image("logo_primerize.png", width='stretch')

st.sidebar.markdown("<br>", unsafe_allow_html=True)
st.sidebar.header("Baseline 1D Parameters")
sb_prefix = st.sidebar.text_input("Construct Prefix/Name", value="my_rna_construct")
sb_min_tm = st.sidebar.slider("Minimum Overlap Tm (°C)", min_value=45.0, max_value=85.0, value=60.0, step=0.5)

sb_limit_primers = st.sidebar.checkbox("Enforce Exact Primer Count Limit", value=False)
sb_num_primers = None # Changed from 0 to None to satisfy the Primerize backend which contradicts the documentation's claim about 0 being the default for no limit. None is the actual default that disables the limit.
if sb_limit_primers:
    sb_num_primers = st.sidebar.number_input("Exact Number of Primers (Even integers only)", min_value=2, max_value=100, value=8, step=2)

sb_min_len = st.sidebar.number_input("Minimum Primer Length (nt)", min_value=10, max_value=50, value=15)
sb_max_len = st.sidebar.number_input("Maximum Primer Length (nt)", min_value=40, max_value=120, value=60)

# Sidebar note about sequence length limits for free tier stability and performance
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

# ---- TAB 1: 1D SIMPLE ASSEMBLY PRIMER DESIGN ----
with tab1:
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
                        
                        # Save execution state inside the browser session storage for child dependencies
                        st.session_state['active_job_1d'] = job_1d
                        
                        # Render Data View Layout Table
                        st.subheader("Primary Assembly Primer Set")
                        primer_records = [{"Primer Index": i + 1, "Oligo Sequence (5' → 3')": s} for i, s in enumerate(job_1d.primer_set)]
                        st.dataframe(pd.DataFrame(primer_records), width="stretch")
                        
                        # Generate the IDT Order Capture Form Output
                        st.subheader("Full Backend Log")
                        st.caption("Raw output engine execution log.")
                        
                        # Safely trap native output print strings
                        buffer = io.StringIO()
                        sys.stdout = buffer
                        print(job_1d)
                        sys.stdout = sys.__stdout__
                        output_text = buffer.getvalue()
                        
                        # Process output text through our clean transformer engine
                        clean_html_output = convert_ansi_to_html(output_text)

                        # Render inside a dedicated CSS flex-grow viewport block
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

# ---- TAB 2: 2D CHEMICAL MAPPING LIBRARY ----
with tab2:
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

                        # Process output text through our clean transformer engine
                        clean_html_output_2d = convert_ansi_to_html(output_text_2d)

                        # Render inside a dedicated CSS flex-grow viewport block
                        st.markdown(
                            f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 18px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d; letter-spacing: 0.5px;">{clean_html_output_2d}</div>', 
                            unsafe_allow_html=True
                        )
                        st.markdown("<br>", unsafe_allow_html=True)
                        
                        col_dl1, col_dl2 = st.columns(2)
                        with col_dl1:
                            st.download_button("Download 2D Plate Specification File (.txt)", data=strip_ansi(output_text_2d), file_name=f"{st.session_state['active_job_1d'].name}_2D_library.txt")
                            
                        with col_dl2:
                            import glob
                            # Generate and serve the Excel file
                            xls_name = f"{st.session_state['active_job_1d'].name}_2D_library"
                            job_2d.save('table', path='.', name=xls_name)
                            xls_files = glob.glob(f"{xls_name}_plate_*.xls")
                            for i, xls_path in enumerate(xls_files):
                                with open(xls_path, "rb") as f:
                                    xls_data = f.read()
                                st.download_button(f"Download Plate {i+1} (.xls)", data=xls_data, file_name=xls_path, mime="application/vnd.ms-excel", key=f"dl_2d_{i}")
                                os.remove(xls_path) # Cleanup
                                
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

# ---- TAB 3: 3D STRUCTURE MUTATIONS ----
with tab3:
    st.header("3D Structure Mutations")
    st.caption("Optimize assembly primer sets for structural targets by integrating thermodynamic calculations with dot-bracket folding models.")

    # 1. ENFORCE GATEKEEPER: Session State Requirement from Tab 1
    if 'active_job_1d' not in st.session_state:
        st.warning(
            "**Prerequisite Required:** You must run a 1D Baseline calculation in Tab 1 "
            "to unlock the 3D structure-guided modules. This populates the foundational "
            "primer landscape required for mutation mapping."
        )
    else:
        # Retrieve the baseline calculation structure from session memory
        job_1d = st.session_state['active_job_1d']
        
        st.info(f"Active Baseline: {job_1d.name} ({len(job_1d.sequence)} nt)")
        
        # 2. SEED THE TARGET SEQUENCE AUTOMATICALLY FOR USER CONVENIENCE
        seq_3d = st.text_area(
            "Target Sequence (Synchronized from Tab 1)", 
            value=job_1d.sequence,
            disabled=True,  # Keeps it locked to match the 1D project exactly
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
        
        # 3. RUN PIPELINE TRIGGER BUTTON
        if st.button("Generate 3D Mutation Design", type="primary", key="btn_3d"):
            st.session_state['run_3d'] = True
            st.session_state['t3_struct'] = struct_3d
            
        if st.session_state.get('run_3d', False):
            cleaned_struct = "".join(st.session_state['t3_struct'].split())
            
            # --- VALIDATION LAYER ---
            if not cleaned_struct:
                st.error("Submission Failure: The Dot-Bracket Structure map cannot be blank.")
                
            # Check that lengths match exactly character-for-character
            elif len(job_1d.sequence) != len(cleaned_struct):
                st.error(
                    f"**Input Length Mismatch!**\n\n"
                    f"The 3D structural folding map must match your target sequence length character-for-character:\n"
                    f"* Your target sequence is **{len(job_1d.sequence)}** nucleotides long.\n"
                    f"* Your dot-bracket structure input is **{len(cleaned_struct)}** characters long.\n\n"
                    f"Please adjust your structure layout to align perfectly before running the pipeline."
                )
            else:
                # --- EXECUTION LAYER ---
                with st.spinner("Calculating 3D structure-guided mutations and primer boundaries..."):
                    try:
                        # Instantiate the 3D module on the singleton worker
                        prm_3d = primerize.Primerize_3D
                        prm_3d.reset()
                        
                        # Set default structural variables to worker environment profile
                        prm_3d.set('prefix', str(job_1d.name))
                        
                        # Execute backend algorithm calculation matrix pass
                        # FIXED: Wrapping cleaned_struct in a list [] and using keyword 'structures'
                        job_3d = prm_3d.design(
                            sequence=job_1d.sequence,
                            primer_set=job_1d.primer_set,
                            structures=[cleaned_struct],  # Documentation expects list(str)
                            prefix=str(job_1d.name)
                        )
                        
                        if job_3d.is_success:
                            st.success("3D Structural Synthesis Optimization Completed.")
                            
                            # Render Data View Layout Table using modernized properties
                            st.subheader("Structural Mutation Assembly Primers")
                            primer_records_3d = [{"Primer Index": i + 1, "Oligo Sequence (5' → 3')": s} for i, s in enumerate(job_3d.primer_set)]
                            st.dataframe(pd.DataFrame(primer_records_3d), width="stretch")
                            
                            # Intercept printed alignment layout from standard streams
                            buffer_3d = io.StringIO()
                            sys.stdout = buffer_3d
                            print(job_3d)
                            sys.stdout = sys.__stdout__
                            output_text_3d = buffer_3d.getvalue()
                            
                            # Convert escape strings to web-native span markup using the updated function
                            clean_html_output_3d = convert_ansi_to_html(output_text_3d)
                            
                            # Render output inside our scroll-locked container wrapper
                            st.subheader("3D Alignment Matrix Mapping")
                            st.markdown(
                                f'<div style="background-color: #1e1e1e; color: #ffffff; padding: 15px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; white-space: pre; overflow-x: auto; font-size: 13px; line-height: 1.5; border: 1px solid #2d2d2d;">{clean_html_output_3d}</div>', 
                                unsafe_allow_html=True
                            )
                            st.markdown("<br>", unsafe_allow_html=True)
                            
                            st.markdown("<br>", unsafe_allow_html=True)
                            
                            col_dl1, col_dl2 = st.columns(2)
                            with col_dl1:
                                st.download_button(
                                    "Download 3D Design Matrix (.txt)", 
                                    data=strip_ansi(output_text_3d), 
                                    file_name=f"{job_1d.name}_3D_structural_design.txt"
                                )
                                
                            with col_dl2:
                                import glob
                                # Generate and serve the Excel file
                                xls_name = f"{job_1d.name}_3D_structural_design"
                                job_3d.save('table', path='.', name=xls_name)
                                xls_files = glob.glob(f"{xls_name}_plate_*.xls")
                                for i, xls_path in enumerate(xls_files):
                                    with open(xls_path, "rb") as f:
                                        xls_data = f.read()
                                    st.download_button(f"Download Plate {i+1} (.xls)", data=xls_data, file_name=xls_path, mime="application/vnd.ms-excel", key=f"dl_3d_{i}")
                                    os.remove(xls_path) # Cleanup
                                    
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
                                    
                        else:
                            st.error("Error: Primerize Engine Failure. No valid primer combinations could satisfy the 3D folding constraints.")
                            
                    except Exception as e:
                        st.error(f"Execution Error Intercepted: {str(e)}")

# ---- TAB 4: CUSTOM CONSTRUCT TARGETED DESIGN ----
with tab4:
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
                    job_custom, error_log = cached_design_custom(st.session_state['active_job_1d'], t4_offset, custom_mut_input)
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
                            import glob
                            # Generate and serve the Excel file
                            xls_name = f"{st.session_state['active_job_1d'].name}_custom"
                            job_custom.save('table', path='.', name=xls_name)
                            xls_files = glob.glob(f"{xls_name}_plate_*.xls")
                            for i, xls_path in enumerate(xls_files):
                                with open(xls_path, "rb") as f:
                                    xls_data = f.read()
                                st.download_button(f"Download Plate {i+1} (.xls)", data=xls_data, file_name=xls_path, mime="application/vnd.ms-excel", key=f"dl_4d_{i}")
                                os.remove(xls_path) # Cleanup
                                
                        with st.expander("Supplementary Data Files"):
                            st.caption("Supplementary output files for downstream analysis and sequence alignment pipelines.")
                            c1, c2, c3 = st.columns(3)
                            constructs_txt, assembly_txt, _ = generate_advanced_files(job_custom)
                            with c1:
                                st.download_button("Download Constructs Map (.txt)", data=constructs_txt, file_name=f"{job_custom.name}_constructs.txt", key="t4_dl_const")
                            with c2:
                                st.download_button("Download Assembly DNA (.txt)", data=assembly_txt, file_name=f"{job_custom.name}_assembly.txt", key="t4_dl_assem")
                                
                    else:
                        clean_err = strip_ansi(error_log).strip()
                        # Extract the specific ValueError message for a cleaner non-technical UX
                        val_err_match = re.search(r'ValueError:\s*(?:ERROR:\s*)?(.*)', clean_err)
                        if val_err_match:
                            parsed_err = val_err_match.group(1).strip()
                            st.error(f"Error: Custom Construction Intercept. The target variations could not map to the active sequence range safely.\n\n**{parsed_err}**")
                        else:
                            st.error(f"Error: Custom Construction Intercept. The target variations could not map to the active sequence range safely.\n\n**Backend Output:**\n`{clean_err}`")
                except Exception as e:
                    st.error(f"Error Compiling Custom Plate Array: {str(e)}")

# ---- TAB 5: PROTOCOLS & ORDERING ----
with tab5:
    st.header("Protocols & Ordering")
    st.markdown("""
For complete and up-to-date documentation regarding the Primerize pipeline, including DNA template design, IDT oligo ordering, PCR assembly protocols, and in vitro transcription (IVT), please refer to the official documentation:

**[Primerize Pipeline Protocols](https://primerize.stanford.edu/protocol/#pipe)**
""")

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