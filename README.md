# Primerize Web Server Suite

This repository provides a modernized, open-source, and portable graphical user interface (GUI) for the **Primerize** platform, originally developed by the Das Lab at Stanford University. 

Primerize utilizes dynamic programming-based algorithms to design optimized primers for PCR assembly, primarily targeting high-throughput RNA synthesis and multidimensional chemical mapping experiments. This interface is built on [Streamlit](https://streamlit.io/) and natively integrates the core `primerize` Python thermodynamics engine, fully replacing the legacy Django/MySQL AWS server architecture decommissioned in May 2026.

## Core Features

- **1D Simple Assembly:** Automated design of PCR assembly primers given a desired DNA/RNA template sequence, minimizing thermodynamic mispriming. Outputs natively formatted text blocks for bulk ordering via IDT DNA Oligos.
- **2D Chemical Mapping:** High-throughput mutagenesis plate array generation for mutate-and-map (M2) multidimensional chemical mapping pipelines.
- **3D Structure Mutations:** Integrates dot-bracket secondary structure folding models to evaluate structural constraints and generate targeted mutation plates (mutate-map-rescue).
- **Custom Construct Mutants:** Allows for explicit targeted variation mapping (e.g., `A10C`) across the active baseline sequence.

## Repository Structure

```text
├── .github/workflows/   # CI/CD deployment configurations (GitHub Actions)
├── .streamlit/          # Web application configuration and UI themes
├── primerize/           # The core Python thermodynamics and algorithm engine
├── app.py               # The main Streamlit web application and UI routing logic
├── Dockerfile           # Docker configuration for isolated environment hosting
├── requirements.txt     # Python package dependencies
└── README.md            # Repository documentation
```

## Local Installation & Usage

To run this application locally, ensure you have Python 3.11+ installed.

1. Clone this repository to your local environment.
2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Launch the Streamlit server:
   ```bash
   streamlit run app.py
   ```
4. Access the web interface at `http://localhost:8501`.

## Official Protocols & Documentation

For complete laboratory protocols regarding DNA template design, IDT oligo ordering, PCR assembly, and in vitro transcription (IVT), please refer to the official documentation:

**[Primerize Pipeline Protocols](https://primerize.stanford.edu/protocol/#pipe)**

## Academic Citation Acknowledgments

**This project is distributed under the MIT license.**

Please cite the use of Primerize with the following foundational publications:

- **Primerize-2D:** Tian, S., and Das, R. (2017) [Automated primer design for RNA multidimensional chemical mapping](https://academic.oup.com/bioinformatics/article-abstract/33/9/1405/2801460/Primerize-2D-automated-primer-design-for-RNA). *Bioinformatics* 33(9): 1405-1406.
- **Primerize:** Tian, S., Yesselman, J.D., Cordero, P., and Das, R. (2015) [Primerize: automated primer assembly for transcribing interesting RNAs](http://nar.oxfordjournals.org/content/43/W1/W522). *Nucleic Acids Research* 43(W1): W522-W526.

### Related Publications on the Mutate-Map-Rescue Pipeline:
- Tian, S., and Das, R. (2016) RNA structure through multidimensional chemical mapping. *Quarterly Review of Biophysics* 49(e7): 1-30.
- Tian, S., Cordero, P., Kladwang, W., and Das, R. (2014) High-throughput mutate-map-rescue evaluates SHAPE-directed RNA structure and uncovers excited states. *RNA* 20(11): 1815-1826.
- Kladwang, W., VanLang, C.C., Cordero P., and Das, R. (2011) A two-dimensional mutate-and-map strategy for non-coding RNA structure. *Nature Chemistry* 3: 954-962.
