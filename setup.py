import setuptools

setuptools.setup(
    name="PEKA",
    version="1.0.0",
    author="Aram Amalietti",
    author_email="aram.amalietti@gmail.com",
    description=(
        "Positionally-enriched k-mer analysis (PEKA) is a software package for "
        "identifying enriched protein-RNA binding motifs from CLIP datasets."
    ),
    long_description=(
        "Positionally-enriched k-mer analysis (PEKA) is a software package for "
        "identifying enriched protein-RNA binding motifs from CLIP datasets. "
        "PEKA compares k-mer enrichment in proximity of high-confidence "
        "crosslink sites (tXn - thresholded crosslinks), located within "
        "crosslinking peaks and having a high cDNA count, relative to "
        "low-count crosslink sites located outside of peaks (oXn - outside "
        "crosslinks). This approach reduces the effects of technical biases, "
        "such as uridine-preference of UV crosslinking. Each k-mer is assigned "
        "a PEKA score, which is used to rank the k-mers from the most to the "
        "least enriched. Additionally, PEKA provides comprehensive "
        "visualisations of motif enrichment profiles around the "
        "high-confidence crosslink sites and clusters the motifs that display "
        "similar profiles. PEKA also enables motif discovery within specific "
        "transcriptomic regions, including or excluding repetitive elements."
    ),
    url="https://github.com/ulelab/peka",
    scripts = [
        "peka.py"
    ],
    python_requires=">=3.7",
    entry_points={"console_scripts": ["peka = peka:main",],},
)
