# Primerize Web 🧬

**Primerize Web** is a web application built on the Primerize Python package for primer design and nucleic acid thermodynamics, developed by the Das Lab at Stanford University for high-throughput RNA synthesis and design.

This website lets you run the Primerize algorithms without writing code. You can design primers for simple 1D assembly, 2D chemical mapping libraries, and 3D structure-guided mutations, and more!

*This website is not affiliated with the official Primerize web server, Stanford University, or the Das Lab.*

## Why build Primerize Web?
The original Primerize web server was officially decommissioned in May 2026. However, the underlying Python algorithmic engine (`primerize`) is still incredibly powerful and widely used in the RNA community.

The code for the original Primerize web server was not released publicly, so we built this repository to resurrect the web server from scratch. The server stack was rebuilt natively in Python using [Streamlit](https://streamlit.io/).

## What is Primerize?
Primerize is a thermodynamic engine that designs optimized PCR assembly primers for high-throughput RNA synthesis. If a researcher wants to map the 2D or 3D structure of an RNA molecule, they need to design dozens or hundreds of targeted mutations. Doing this by hand is a nightmare. Primerize automates it by calculating thermodynamic folding models to minimize mispriming.

## How it Works (Under the Hood)
If you're looking to contribute, here's a primer (pun intended) on how the app is architected:

- **Streamlit:** The entire UI is written in pure Python using Streamlit. This makes it incredibly easy to spin up natively or containerize if you prefer to use docker. Guides detailing how to dockerize a Streamlit app are easy to find online.
- **Modular Tabs:** The UI is broken down into `tabs/` (e.g., `tab1_baseline.py`, `tab2_mapping.py`). This makes it easy to add new features to a specific page without breaking existing ones.
- **RAM-only processes** The only feature that writes to disk is the .xls plate exports on tabs 2, 3, and 4. The other features work entirely in-memory, which helps keep Primerize Web extremely fast.
- **Thread-Safe Concurrency:** The core primerize math engine was originally designed for a single user, meaning that if 50 people visit the website at once, they are all physically sharing the same engine in the background. If two researchers clicked "Primerize!" at the exact same time, their sequences would collide inside the engine and corrupt the results. To make this safe for a public web application, we implemented a queuing lock (`threading.Lock()`). When you submit a sequence, the app temporarily "locks" the engine, does your math, saves your result, and then unlocks it for the next person in line.
- **Temporary Sessions & State:** There is no backend database. User sessions are managed purely through Streamlit's ephemeral `st.session_state` dictionary. Every visitor gets their own isolated state memory tied to their active browser tab; if the tab is closed or refreshed, the server dumps the memory.
- **Cross-Tab Synchronization:** The application enforces a linear workflow. The user must first generate a baseline sequence in Tab 1, which the app saves into `st.session_state['active_job_1d']`. Tabs 2, 3, and 4 read directly from this session state, allowing the user to seamlessly generate complex mapping libraries without ever having to re-enter their sequence.
- **ANSI-to-HTML Log Intercept:** The legacy `primerize` engine outputs 3D alignment matrices directly to the console (`stdout`) using ANSI color codes. Since Python's print() statements output to the server's terminal rather than the user's web page, we built a wrapper in `utils.py` using `contextlib.redirect_stdout` to intercept the print stream, parse the ANSI color codes into raw HTML/CSS, and render the terminal output directly inside the Streamlit UI.

## How to Work on It
We would love your help! All contributions are welcome.

### Getting Started
1. Clone the repository:
   ```bash
   git clone https://github.com/your-org/Primerize_GUI.git
   cd Primerize_GUI
   ```
2. Create a virtual environment and install the dependencies:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows use: .venv\Scripts\activate
   pip install -r requirements.txt
   ```
3. Launch the development server:
   ```bash
   streamlit run app.py
   ```
4. Access the web interface at `http://localhost:8501`.

### Contributing
If you have an idea or see a bug, please feel free to:
1. **Fork the repo** and create your branch (`git checkout -b feature/AmazingFeature`)
2. **Commit your changes** (`git commit -m 'Add some AmazingFeature'`)
3. **Push to the branch** (`git push origin feature/AmazingFeature`)
4. **Open a Pull Request!** We review PRs regularly and are always happy to collaborate and help you get your code merged.

---

## Privacy Policy

**Primerize Web does not collect or track your data.** Primerize Web does not use a database or log any personally identifiable information. We do not log or permenantly store your DNA/RNA sequence data, project names, or generated primers. Session data is strictly tied to your active browser tab and destroyed if you close or reload the tab.

**If you visit [primerize.streamlit.app](primerize.streamlit.app):** Because this application is hosted on the free Streamlit Community Cloud platform, Streamlit (and its parent company, Snowflake) automatically collects some broad usage information. This is standard for cloud hosting providers and generally includes basic server access logs (such as IP addresses and browser types) and anonymous page-view analytics used to maintain their platform security and stability. You can review their full platform privacy policy [here](https://streamlit.io/privacy-policy).

Primerize Web is open-source software under the MIT License. The original Primerize Python package is also released under the MIT License. You can review the code yourself to verify that it does what we say it does.

Organizations with strict internal infosec policies prohibiting the transmission of proprietary sequences to external cloud servers can easily deploy the open-source repository on an isolated internal network.

---

## Official Protocols & Documentation
For complete laboratory protocols regarding DNA template design, IDT oligo ordering, PCR assembly, and in vitro transcription (IVT), please refer to the official documentation:
**[Primerize Pipeline Protocols](https://primerize.stanford.edu/protocol/#pipe)**

## Academic Citation Acknowledgments
**This project is distributed under the MIT license.**

If you use Primerize in your research, please cite the foundational publications:
- **Primerize-2D:** Tian, S., and Das, R. (2017) [Automated primer design for RNA multidimensional chemical mapping](https://academic.oup.com/bioinformatics/article-abstract/33/9/1405/2801460/Primerize-2D-automated-primer-design-for-RNA). *Bioinformatics* 33(9): 1405-1406.
- **Primerize:** Tian, S., Yesselman, J.D., Cordero, P., and Das, R. (2015) [Primerize: automated primer assembly for transcribing interesting RNAs](http://nar.oxfordjournals.org/content/43/W1/W522). *Nucleic Acids Research* 43(W1): W522-W526.

**Related Publications on the Mutate-Map-Rescue Pipeline:**
- Tian, S., and Das, R. (2016) [RNA structure through multidimensional chemical mapping.](http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=10242118&fulltextType=RV&fileId=S0033583516000020) *Quarterly Review of Biophysics* 49(e7): 1-30.
- Tian, S., Cordero, P., Kladwang, W., and Das, R. (2014) [High-throughput mutate-map-rescue evaluates SHAPE-directed RNA structure and uncovers excited states.](http://rnajournal.cshlp.org/content/20/11/1815) *RNA* 20(11): 1815-1826.
- Kladwang, W., VanLang, C.C., Cordero P., and Das, R. (2011) [A two-dimensional mutate-and-map strategy for non-coding RNA structure.](http://www.nature.com/nchem/journal/v3/n12/full/nchem.1176.html) *Nature Chemistry* 3: 954-962.
