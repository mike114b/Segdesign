# SegDesign: Intelligent Protein Segment Design Pipeline

<div align="center">

**An integrated pipeline for intelligent protein segment design combining sequence analysis, structure prediction, and generative modeling**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

## 📖 Overview

SegDesign is an automated pipeline for intelligent protein segment design. It integrates multiple state-of-the-art bioinformatics tools and deep learning models to perform comprehensive protein analysis and design:

- **Sequence Conservation Analysis**: Using HMMER for evolutionary conservation analysis
- **Structure Generation**: Using RFdiffusion for targeted protein backbone generation
- **Sequence Design**: Using ProteinMPNN for amino acid sequence optimization
- **Structure Validation**: Using ESMFold for predicted structure quality assessment
- **Sequence Clustering**: Using MMSeqs2 for sequence similarity analysis

## 🏗️ Architecture

```
SegDesign/
├── Segdesign.py              # Main entry point
├── Segdesign/
│   ├── hmmer/               # Sequence conservation analysis
│   ├── rfdiffusion/         # Structure generation
│   ├── mpnn/                # Sequence design
│   ├── esmfold/             # Structure prediction
│   └── dssp/                # Secondary structure analysis
├── config/
│   ├── config.yaml          # User configuration
│   └── setting.yaml         # System settings
├── environments/            # Environment installation scripts
└── example/                 # Example outputs
```

## 🚀 Quick Start

### Prerequisites

- **Operating System**: Linux (recommended) or Windows with WSL2
- **Python**: 3.9 or higher
- **Conda/Miniconda**: Required for environment management
- **GPU**: NVIDIA GPU with CUDA support (strongly recommended for ESMFold and RFdiffusion)
- **Memory**: At least 16GB RAM (32GB+ recommended)
- **Storage**: At least 200GB free space

### Installation

#### 1. Clone the Repository

```bash
git clone https://github.com/mike114b/SegDesign2.git
cd SegDesign2
```

#### 2. Install Conda Environments

The project requires 3 conda environments to run different modules:
- **segdesign**: Main environment containing HMMER, MMSeqs2, DSSP, etc.
- **segdesign_esmfold**: For running ESMFold model
- **segdesign_SE3nv**: For running RFdiffusion model and ProteinMPNN model

We provide installation scripts in the `environments/` directory for user convenience. Before running the scripts, ensure that Conda or Miniconda is installed. You can use CONDA_PATH to specify the Conda installation path. If not specified, ensure conda runs properly and the script will use the conda run command to install environments.

Please run the installation scripts:

```bash
# You can set CONDA_PATH to specify Conda installation path, but this is optional
CONDA_PATH="/path/to/your/anaconda3"

# Install main environment (HMMER, MMSeqs2, DSSP, etc.)
bash ./environments/segdesign_env.sh

# Install SE3nv environment (containing RFdiffusion and ProteinMPNN)
bash ./environments/segdesign_SE3nv_env.sh

# Install ESMFold environment (requires CUDA support)
bash ./environments/esmfold_env.sh
```

#### 3. Install Databases (Optional)

For HMMER analysis, you may need to download sequence databases:

```bash
# Download UniRef90 database
bash environments/download_uniref90.sh

# Download UniRef100 database
bash environments/download_uniref100.sh
```

#### 4. Configure Paths

You can edit `config/setting.yaml` file to configure the following paths:
- Anaconda installation path
- RFdiffusion installation path
- ProteinMPNN installation path
In general, you don't need to modify these paths, using default values is fine.

## 📋 Configuration

### User Configuration (`config/config.yaml`)

The user configuration file controls the workflow parameters:

```yaml
project:
  anaconda_path:                     # Anaconda installation path, leave empty to use conda run command
  input_pdb: ./Dusp4.pdb             # Input protein structure file
  output_dir: ./output               # Output directory
  chain: A                           # Chain to analyze
  sequence_length: 394               # Full sequence length
  segment: 346-394                   # Design region (optional)

profile:
  database: ./uniprot_sprot.fasta    # Sequence database
  bitscore: 0.3                      # HMMER bit score threshold
  n_iter: 5                          # JackHMMER iterations
  cpu: 10                            # Number of CPU cores
  threshold: 0.6                     # Conservation threshold

rfdiffusion:
  num_designs: 10                    # Number of designs to generate
  threshold: 0.04                    # Design quality threshold
  helix: false                       # Design as alpha-helix
  strand: false                      # Design as beta-strand

mpnn:
  num_seq_per_target: 20             # Sequences per design
  sampling_temp: 0.3                 # MPNN sampling temperature
  seed: 42                           # Random seed
  top_percent: 0.9                   # Top percentage selection

esmfold:
  ptm_threshold: 0.54                # PTM score threshold
  plddt_threshold: 70                # pLDDT score threshold
```

## 💻 Usage

### Basic Usage

Run the complete pipeline:

```bash
python Segdesign.py --config config/config.yaml
```

### Individual Module Execution

Individual modules can be run separately:

```bash
# Run sequence analysis only
 conda run -n segdesign python ./SegDesign2/Segdesign/hmmer/hmmer.py --input_pdb ./Dusp4.pdb --select_chain A --output_folder ./Dusp4_example/hmmer_out --bitscore 0.3 --n_iter 5 --database ./uniprot_sprot.fasta --cpu 10 --minimum_sequence_coverage 50 --minimum_column_coverage 70 --final_report_folder ./Dusp4_example


# Run protein backbone design only
conda run -n segdesign_SE3nv python /home/wangxuming/SegDesign2_test/Segdesign/rfdiffusion/rf_diffusion.py --run_inference_path ./RFdiffusion/scripts/run_inference.py --inference.input_pdb ./Dusp4.pdb --inference.output_prefix ./Dusp4_example/rfdiffusion_out/sample/Dusp4_A --inference.num_designs 10 --contigmap.contigs '[A1-394]' --contigmap.inpaint_str '[A346-394]' --diffuser.partial_T 50 --contigmap.inpaint_str_strand '[A346-394]'

# Run structure prediction only
conda run -n segdesign_esmfold python ./SegDesign2/Segdesign/esmfold/esmfold.py --input_pdb ./Dusp4.pdb --output_folder ./Dusp4_example/esmfold_out --ptm_threshold 0.54 --plddt_threshold 70
```

### Example: Dusp4 Protein Design

The `example/Dusp4_example/` directory contains a complete output example:

```bash
# Run the example workflow
python Segdesign.py --config example/Dusp4_example/config.yaml
```

## 📊 Output Structure

```
output/
├── config.yaml                    # Copy of configuration
├── hmmer_out/                     # HMMER analysis results
│   ├── Dusp4_A_Recommended_Design_Area.txt
│   ├── Dusp4_A_conservative_comprehensive_report.csv
│   └── jackhmmer_out/            # Raw HMMER alignments
├── rfdiffusion_out/              # RFdiffusion results
│   ├── sample/                   # Generated backbones
│   └── filter_results/           # Filtered structures
├── mpnn_out/                     # MPNN sequence designs
│   ├── seqs/                     # Designed sequences
│   └── csv_files/                # Analysis CSVs
└── esmfold_report.csv            # Final validation report
```

### Output Columns Description

| Column | Description |
|--------|-------------|
| index | Design identifier |
| backbone | Source backbone structure |
| segment | Designed region |
| score | Design score |
| plddt_score | ESMFold pLDDT confidence score |
| ptm_score | ESMFold PTM score |
| whether_pass | Quality control pass status |

## 🔧 Module Details

### 1. HMMER Module
- Performs sequence conservation analysis using JackHMMER
- Identifies conserved regions for intelligent design area selection
- Generates comprehensive conservation reports

### 2. RFdiffusion Module
- Generates novel protein backbones for the design region
- Supports secondary structure constraints (helix/strand)
- Produces multiple design candidates

### 3. ProteinMPNN Module
- Designs amino acid sequences for generated backbones
- Optimizes sequences for stability and expression
- Supports fixed backbone positions

### 4. ESMFold Module
- Validates designed structures using deep learning prediction
- Assesses pLDDT and PTM scores
- Filters low-quality designs

### 5. MMSeqs2 Module (Optional)
- Performs sequence clustering analysis
- Identifies sequence diversity
- Generates cluster reports

## ⚠️ Troubleshooting

### GPU Memory Issues
```bash
# Reduce batch size or number of designs
# Set environment variable for GPU memory limit
export CUDA_VISIBLE_DEVICES=0
```

### Conda Environment Activation
```bash
# Ensure CONDA_PATH is set correctly
export CONDA_PATH="/path/to/anaconda3"
source $CONDA_PATH/etc/profile.d/conda.sh
```

### Database Errors
- Verify database paths in `config/setting.yaml`
- Ensure databases are properly formatted
- Check file permissions

## 📧 Contact

For questions or suggestions, please open an issue or contact the author.

---

<div align="center">

**Happy Protein Designing! 🔬🧬**

</div>
