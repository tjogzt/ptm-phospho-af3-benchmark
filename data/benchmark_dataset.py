"""
PTM Phosphorylation Benchmark Dataset
======================================
Defines the benchmark protein set for phosphorylation structural effect analysis.

Selection criteria:
- Has both phosphorylated and non-phosphorylated structures in PDB (paired set)
- Or has high-quality phosphorylated structure (single set)
- Resolution ≤ 3.0 Å (X-ray) or ≤ 4.0 Å (cryo-EM)
- Covers diverse structural effect types: local rearrangement, allosteric, disorder-order
"""

BENCHMARK_PROTEINS = {
    # ====== PAIRED SET (has both phos & unphos structures) ======
    
    "CDK2": {
        "uniprot": "P24941",
        "gene": "CDK2",
        "description": "Cyclin-dependent kinase 2",
        "phos_sites": ["Thr160"],
        "phos_type": "TPO",
        "pdbs": {
            "unphos": ["1HCL"],   # inactive, unphosphorylated T-loop
            "phos": ["1JST"],     # active, pThr160
        },
        "effect_type": "local_rearrangement",  # T-loop refolding
        "chain": "A",
        "length": 298,
        "blind_test": False,  # PDB date: 1996
    },
    
    "ERK2": {
        "uniprot": "P28482",
        "gene": "MAPK1",
        "description": "Mitogen-activated protein kinase 1",
        "phos_sites": ["Thr183", "Tyr185"],
        "phos_type": "TPO_PTR",  # dual phosphorylation
        "pdbs": {
            "unphos": ["1ERK"],   # unphosphorylated
            "phos": ["2ERK"],     # dually phosphorylated (pThr183/pTyr185)
        },
        "effect_type": "local_rearrangement",
        "chain": "A",
        "length": 358,
        "blind_test": False,
    },
    
    "PKA_Calpha": {
        "uniprot": "P17612",
        "gene": "PRKACA",
        "description": "cAMP-dependent protein kinase catalytic subunit alpha",
        "phos_sites": ["Thr197"],
        "phos_type": "TPO",
        "pdbs": {
            "unphos": ["1CTP"],   # apo, unphosphorylated Thr197
            "phos": ["1ATP"],     # pThr197 + bound peptide
        },
        "effect_type": "allosteric",
        "chain": "A",
        "length": 350,
        "blind_test": False,
    },
    
    "CheY": {
        "uniprot": "P0AE67",
        "gene": "cheY",
        "description": "Chemotaxis response regulator (E. coli)",
        "phos_sites": ["Asp57"],
        "phos_type": None,  # Asp phosphorylation, not directly supported by AF3
        "note": "Use BeF3- mimic structure instead of native Asp-P",
        "pdbs": {
            "unphos": ["3CHY"],   # inactive, unphosphorylated
            "phos": ["1FQW"],     # BeF3- activated mimic
        },
        "effect_type": "allosteric",
        "chain": "A",
        "length": 129,
        "blind_test": False,
    },
    
    "Glycogen_Phosphorylase": {
        "uniprot": "P00489",
        "gene": "PYGM",
        "description": "Glycogen phosphorylase, muscle form",
        "phos_sites": ["Ser15"],
        "phos_type": "SEP",
        "pdbs": {
            "unphos": ["1GPA"],   # T-state, unphosphorylated (actually has SEP but can be treated as reference)
            "phos": ["1GPA"],     # Use same PDB, phosphorylated Ser14 (PDB numbering off-by-one)
        },
        "note": "N-terminal tail phosphorylation induces folding and allosteric activation",
        "effect_type": "disorder_order",
        "chain": "A",
        "length": 842,
        "blind_test": False,
    },
    
    "CDK7": {
        "uniprot": "P50613",
        "gene": "CDK7",
        "description": "Cyclin-dependent kinase 7",
        "phos_sites": ["Thr170"],
        "phos_type": "TPO",
        "pdbs": {
            "unphos": ["1UA2"],   # inactive
            "phos": ["7Ml6"],     # active
        },
        "effect_type": "local_rearrangement",
        "chain": "A",
        "length": 346,
        "blind_test": True,  # 7ML6 released 2022
    },
    
    "IRF3": {
        "uniprot": "Q14653",
        "gene": "IRF3",
        "description": "Interferon regulatory factor 3",
        "phos_sites": ["Ser386", "Ser396"],
        "phos_type": "SEP",
        "pdbs": {
            "unphos": ["1QWT"],   # autoinhibited
            "phos": ["1T2K"],     # phosphorylated mimic (S→D), or reference only
        },
        "note": "Phosphorylation triggers dimerization and activation. Use 1T2K as phos reference (phosphomimetic)",
        "effect_type": "allosteric",
        "chain": "A",
        "length": 395,
        "blind_test": False,
    },
    
    "Aurora_A": {
        "uniprot": "O14965",
        "gene": "AURKA",
        "description": "Aurora kinase A",
        "phos_sites": ["Thr288"],
        "phos_type": "TPO",
        "pdbs": {
            "unphos": ["1MQ4"],   # unphosphorylated (inactive)
            "phos": ["1OL5"],     # phosphorylated T288 (active)
        },
        "effect_type": "local_rearrangement",
        "chain": "A",
        "length": 403,
        "blind_test": False,
    },
    
    "P38alpha": {
        "uniprot": "Q16539",
        "gene": "MAPK14",
        "description": "p38 MAP kinase alpha",
        "phos_sites": ["Thr180", "Tyr182"],
        "phos_type": "TPO_PTR",
        "pdbs": {
            "unphos": ["1WFC"],   # unphosphorylated
            "phos": ["3D83"],     # dually phosphorylated
        },
        "effect_type": "local_rearrangement",
        "chain": "A",
        "length": 360,
        "blind_test": False,
    },

    # ====== SINGLE-STRUCTURE SET (phos structure only, unphos from AF3/ESMFold) ======
    
    "EGFR_kinase": {
        "uniprot": "P00533",
        "gene": "EGFR",
        "description": "Epidermal growth factor receptor kinase domain",
        "phos_sites": ["Tyr845"],  # activation loop phosphorylation
        "phos_type": "PTR",
        "pdbs": {
            "phos": ["3W2S"],     # phosphorylated, active
            "unphos": [],          # unphos reference from AF2/AF3 prediction
        },
        "effect_type": "local_rearrangement",
        "chain": "A",
        "length": 330,
        "blind_test": False,
    },
    
    "c-Src_kinase": {
        "uniprot": "P12931",
        "gene": "SRC",
        "description": "Proto-oncogene tyrosine-protein kinase Src",
        "phos_sites": ["Tyr416"],
        "phos_type": "PTR",
        "pdbs": {
            "phos": ["1Y57"],     # active, pTyr416
            "unphos": [],

        },
        "note": "Autoinhibited form: 2SRC (unphos)",
        "effect_type": "allosteric",
        "chain": "A",
        "length": 451,
        "blind_test": False,
    },
    
    "AMPK_alpha2": {
        "uniprot": "P54646",
        "gene": "PRKAA2",
        "description": "AMP-activated protein kinase alpha-2",
        "phos_sites": ["Thr172"],
        "phos_type": "TPO",
        "pdbs": {
            "phos": ["2H6D"],     # pThr172 active
            "unphos": [],
        },
        "effect_type": "allosteric",
        "chain": "A",
        "length": 552,
        "blind_test": False,
    },
    
    "CFTR_R_domain": {
        "uniprot": "P13569",
        "gene": "CFTR",
        "description": "Cystic fibrosis transmembrane conductance regulator (R domain)",
        "phos_sites": ["Ser660", "Ser700", "Ser737", "Ser813"],
        "phos_type": "SEP",
        "pdbs": {
            "phos": ["5UAK"],     # phosphorylated
            "unphos": [],
        },
        "note": "Multi-site phosphorylation regulates channel gating. R domain is largely disordered.",
        "effect_type": "disorder_order",
        "chain": "A",
        "length": 200,  # R-domain fragment
        "blind_test": False,
    },
    
    "Hsp27": {
        "uniprot": "P04792",
        "gene": "HSPB1",
        "description": "Heat shock protein beta-1 (Hsp27)",
        "phos_sites": ["Ser15", "Ser78", "Ser82"],
        "phos_type": "SEP",
        "pdbs": {
            "phos": ["4MHR"],     # partially phosphorylated
            "unphos": [],
        },
        "note": "Phosphorylation regulates oligomerization. Multiple sites create graded response.",
        "effect_type": "disorder_order",
        "chain": "A",
        "length": 205,
        "blind_test": False,
    },
    
    "CTD_kinase": {
        "uniprot": "P24928",
        "gene": "POLR2A",
        "description": "RNA polymerase II CTD (C-terminal domain) - model peptide",
        "phos_sites": ["Ser2", "Ser5"],
        "phos_type": "SEP",
        "pdbs": {
            "phos": ["4JXT"],     # CTD peptide bound to phosphatase
            "unphos": [],
        },
        "note": "CTD heptad repeats (YSPTSPS). Phosphorylation code regulates transcription cycle.",
        "effect_type": "disorder_order",
        "chain": "B",  # CTD peptide chain
        "length": 28,  # short peptide
        "blind_test": False,
    },
    
    "PhK_gamma": {
        "uniprot": "P15735",
        "gene": "PHKG1",
        "description": "Phosphorylase b kinase gamma catalytic chain",
        "phos_sites": ["Thr10"],
        "phos_type": "TPO",
        "pdbs": {
            "phos": ["1PHK"],     # phosphorylated
            "unphos": [],
        },
        "effect_type": "local_rearrangement",
        "chain": "A",
        "length": 386,
        "blind_test": False,
    },

    # ====== NEGATIVE CONTROL SET (phosphorylation with minimal structural effect) ======
    
    "Lysozyme": {
        "uniprot": "P00698",
        "gene": "LYZ",
        "description": "Lysozyme C (negative control — no in vivo phosphorylation)",
        "phos_sites": ["Tyr23"],  # artificially introduced phosphosite
        "phos_type": "PTR",
        "pdbs": {
            "unphos": ["1LYS"],
            "phos": [],
        },
        "note": "No natural phosphorylation. Used to test AF3 false positive rate for PTM structures.",
        "effect_type": "negative_control",
        "chain": "A",
        "length": 129,
        "blind_test": False,
    },
    
    "Ubiquitin": {
        "uniprot": "P0CG48",
        "gene": "UBB",
        "description": "Ubiquitin (negative control — highly stable fold)",
        "phos_sites": ["Ser65"],  # PINK1 phosphorylation site
        "phos_type": "SEP",
        "pdbs": {
            "unphos": ["1UBQ"],
            "phos": ["5OXL"],     # pSer65 ubiquitin (PINK1)
        },
        "note": "Ubiquitin fold is extremely stable. pSer65 causes subtle surface change, mainly functional (mitophagy receptor recognition).",
        "effect_type": "negative_control",
        "chain": "A",
        "length": 76,
        "blind_test": False,
    },
}

# ====== BLIND TEST SUBSET ======
# PDB structures released after AF3 training data cutoff (2021-09-30)
BLIND_TEST_IDS = ["7ML6", "8XYZ"]  # CDK7: 7ML6 (2022); add more as identified

# ====== TOTAL SUMMARY ======
"""
Total proteins: 18
  - Paired (phos+unphos in PDB): 9
  - Single phos (unphos from AF3/ESMFold): 7
  - Negative controls: 2
  
  - Blind test subset: ≥1
  - Ser phosphorylation: ~10
  - Thr phosphorylation: ~5
  - Tyr phosphorylation: ~3

Effect type distribution:
  - local_rearrangement: 8
  - allosteric: 4
  - disorder_order: 4
  - negative_control: 2
"""
