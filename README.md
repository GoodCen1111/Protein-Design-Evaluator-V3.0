# ProteinDesignEvaluator

Advanced Protein Design Comprehensive Evaluation Tool / 蛋白质设计综合评估工具

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8+-green.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

## Overview / 项目简介

ProteinDesignEvaluator is a professional evaluation software for protein design research, integrating comprehensive analysis of protein sequences and structures. The tool provides detailed physicochemical property analysis, post-translational modification (PTM) site prediction, and structural quality assessment.

蛋白质设计综合评估工具是一款面向蛋白质设计研究领域的专业评估软件，集成了蛋白质序列和结构的综合分析功能。

## Features / 主要功能

### Sequence Analysis / 序列分析
- **Isoelectric Point (pI)** - Affects purification, electrophoresis, solubility
- **Hydrophobicity Profile** - Detect strong hydrophobic patches
- **Antigenicity Prediction** - Kolaskar & Tongaonkar method
- **Flexibility/Rigidity Distribution** - Karplus & Schulz method
- **Instability Index** - Guruprasad method for protein stability
- **Aliphatic Index** - Thermal stability assessment
- **Transmembrane Region Prediction** - Kyte-Doolittle method

### Post-Translational Modifications / 翻译后修饰
- **Disulfide Bond Prediction** - Cysteine pairing analysis
- **N-Glycosylation Sites** - Asn-X-Ser/Thr motif
- **O-Glycosylation Sites** - Ser/Thr residue analysis
- **Phosphorylation Sites** - S/T/Y phosphorylation potential

### Structural Analysis / 结构分析
- **PDB Structure Analysis** - Multi-chain automatic recognition
- **Geometric Parameter Calculation** - Peptide bonds, CA steps, dihedral angles
- **Clash Detection** - Spatial collision detection
- **Ramachandran Plot** - Backbone conformation visualization

### Visualization / 可视化
- Hydrophobicity distribution plots
- Flexibility/Rigidity profiles
- Antigenicity profiles
- PTM site maps
- Comprehensive analysis dashboard
- Batch comparison charts

## Installation / 安装

### Prerequisites / 前置要求
- Python 3.8 or higher
- Windows 7/10/11, Linux, or macOS

### Option 1: From Source / 方式一：从源码运行

```bash
# Clone or download the repository
cd ProteinDesignEvaluator

# Install dependencies
pip install -r requirements.txt

# Run the application
python main.py
```

### Option 2: Pre-built Executable / 方式二：使用预编译程序

Download the appropriate executable for your platform from the releases page:

- **Windows**: `ProteinDesignEvaluator.exe`
- **Linux**: `ProteinDesignEvaluator` (no installation required)

> Note: The executable version requires no Python installation.

## Usage / 使用说明

### Sequence Analysis / 序列分析

1. Navigate to the "Sequence Analysis" tab
2. Enter amino acid sequence (e.g., `MVLSPADKTN...`)
3. Click "Analyze" to get comprehensive results

### PDB Structure Analysis / PDB结构分析

1. Navigate to the "PDB Analysis" tab
2. Select a PDB file or directory containing PDB files
3. Click "Analyze" to evaluate structure quality

### Comprehensive Analysis / 综合分析

1. Provide both PDB and sequence data
2. Generate comprehensive reports and visualizations

## Output / 输出

Analysis results are automatically saved to:
- **Windows**: `~/ProteinDesignEvaluator_Results/`
- **Linux/macOS**: `~/ProteinDesignEvaluator_Results/`

Output formats:
- PNG (charts and plots)
- CSV (tabular data)
- JSON (complete analysis data)
- TXT (human-readable reports)

## Scoring System / 评分系统

| Score Range | Grade | Description |
|------------|-------|-------------|
| 80-100 | Excellent | High quality, suitable for experimental validation |
| 60-79 | Good | Good quality, can be optimized |
| 40-59 | Average | Some issues, careful selection needed |
| 0-39 | Poor | Significant problems, redesign recommended |

## Project Structure / 项目结构

```
ProteinDesignEvaluator/
├── main.py                 # Main application entry
├── src/
│   ├── __init__.py
│   ├── analyzer.py         # Sequence analysis module
│   ├── pdb_analyzer.py     # PDB structure analysis
│   └── visualizer.py       # Visualization module
├── requirements.txt         # Python dependencies
├── ProteinDesignEvaluator.spec  # PyInstaller configuration
├── build_exe.py            # Build script
├── README.md               # This file
└── docs/
    └── USER_GUIDE.md       # Detailed user guide
```

## Dependencies / 依赖

- **PyQt5** >= 5.15.0 - GUI framework
- **numpy** >= 1.20.0 - Numerical computation
- **matplotlib** >= 3.3.0 - Visualization

## Building from Source / 从源码编译

### Windows

```bash
pip install -r requirements.txt
pip install pyinstaller
python build_exe.py
```

The executable will be in `dist/ProteinDesignEvaluator/`

### Linux

```bash
pip install -r requirements.txt
pip install pyinstaller
python build_exe.py
```

## Contributing / 贡献

Contributions are welcome! Please feel free to submit issues or pull requests.

## License / 许可证

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments / 致谢

- Kyte-Doolittle hydrophobicity scale
- Kolaskar & Tongaonkar antigenicity prediction
- Guruprasad instability index method
- Karplus & Schulz flexibility prediction

## Citation / 引用

If you use this tool in your research, please cite:

```
ProteinDesignEvaluator v3.0
Advanced Protein Design Comprehensive Evaluation Tool
```

---

**Version**: 3.0.0
**Last Updated**: 2026-06-11
