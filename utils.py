import streamlit as st
import threading
import io
import re
import sys
import zipfile
from contextlib import redirect_stdout
from ansi2html import Ansi2HTMLConverter
import primerize

# Pre-configure matplotlib for headless server-side rendering
import matplotlib
matplotlib.use('Agg')

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

# Cache execution results to prevent redundant calculations across UI re-renders
@st.cache_data(max_entries=128, ttl=3600)
def cached_design_1d(seq, min_tm, num_primers, min_len, max_len, prefix):
    prm_1d = primerize.Primerize_1D
    
    # 1. Reset the singleton back to factory defaults to clear out old runs
    prm_1d.reset()
    
    # 2. Pass the parameters directly into the design function call execution layer.
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
        job = primerize.Primerize_2D.design(_job_1d, offset=int(offset), which_muts=mut_range, which_lib=int(which_lib))
    return job, *generate_plate_svgs(job)

@st.cache_data(max_entries=64, ttl=1800)
def cached_design_3d(_job_1d, offset, structures, n_mutations, which_lib, is_single, is_fillwt):
    with primerize_lock:
        job = primerize.Primerize_3D.design(
            _job_1d, 
            offset=int(offset), 
            structures=structures, 
            N_mutations=int(n_mutations), 
            which_lib=int(which_lib), 
            is_single=is_single, 
            is_fillWT=is_fillwt
        )
    return job, *generate_plate_svgs(job)

@st.cache_data(max_entries=64, ttl=1800)
def cached_design_custom(_job_1d, offset, raw_mutation_string):
    mut_list = primerize.Construct_List()
    parsed_muts = [m.strip() for m in re.split(r'[,\s;]+', raw_mutation_string) if m.strip()]
    mut_list.push(parsed_muts)
    
    f = io.StringIO()
    with primerize_lock:
        with redirect_stdout(f):
            job = primerize.Primerize_Custom.design(_job_1d, offset=int(offset), mut_list=mut_list)
            
    return job, f.getvalue(), *generate_plate_svgs(job)

def generate_plate_svgs(job, ref_primer=''):
    """Generates SVGs of 96-well plate layouts and zips them in memory."""
    plates_data = [] 
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        for p_idx, primer_row in enumerate(job._data.get('plates', [])):
            for pl_idx, plate in enumerate(primer_row):
                if not plate: continue # Skip empty plates
                
                suffix = ' R' if p_idx % 2 else ' F'
                p_name, pl_name = f"Primer {p_idx + 1}{suffix}", f"Plate {pl_idx + 1}"
                fname = f"Primer_{p_idx + 1}{suffix.replace(' ', '')}_Plate_{pl_idx + 1}.svg"
                
                # Generate SVG in memory
                buf = io.BytesIO()
                plate.save(ref_primer=ref_primer, file_name=buf, title=f"{p_name} | {pl_name}")
                
                # Instantly strip hardcoded dimensions for responsive Streamlit rendering
                svg = re.sub(r'(width|height)="[^"]+"', r'\1="100%"', buf.getvalue().decode('utf-8'), count=2)
                
                plates_data.append({"primer_name": p_name, "plate_name": pl_name, "svg_string": svg})
                zf.writestr(fname, svg)
                    
    return plates_data, zip_buffer.getvalue()
