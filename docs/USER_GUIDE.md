# ProteinDesignEvaluator v3.0 - User Guide / 使用说明

## Table of Contents / 目录

1. [Quick Start / 快速开始](#quick-start--快速开始)
2. [Interface Overview / 界面概览](#interface-overview--界面概览)
3. [Sequence Analysis / 序列分析](#sequence-analysis--序列分析)
4. [PDB Structure Analysis / PDB结构分析](#pdb-structure-analysis--pdb结构分析)
5. [Comprehensive Analysis / 综合分析](#comprehensive-analysis--综合分析)
6. [Results Export / 结果导出](#results-export--结果导出)
7. [FAQ / 常见问题](#faq--常见问题)

---

## Quick Start / 快速开始

### Windows Executable / Windows可执行程序

1. Download `ProteinDesignEvaluator.exe`
2. Double-click to run (no installation required)
3. Start analyzing!

### From Source / 从源码运行

```bash
pip install -r requirements.txt
python main.py
```

---

## Interface Overview / 界面概览

The application has 5 main tabs:

| Tab | Function |
|-----|----------|
| Sequence Analysis | Analyze amino acid sequences |
| PDB Structure Analysis | Analyze protein structure files |
| Comprehensive Analysis | Combined sequence + structure analysis |
| Results View | View and export results |
| Help | Usage instructions |

---

## Sequence Analysis / 序列分析

### Single Sequence Analysis / 单序列分析

1. Enter sequence name (optional)
2. Enter amino acid sequence
3. Click **"Analyze Sequence"**

### Supported Amino Acids / 支持的氨基酸

Standard 20 amino acids: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V

### Batch Analysis / 批量分析

1. Click **"Select FASTA"** to choose a FASTA file
2. Click **"Batch Analysis"**

### Analysis Results / 分析结果

The tool provides:

- **Overall Score** (0-100) with grade (Excellent/Good/Average/Poor)
- **Isoelectric Point (pI)** - affects purification and solubility
- **Molecular Weight** - in kDa
- **Hydrophobicity** - average and profile
- **Instability Index** - II>40 indicates unstable protein
- **Aliphatic Index** - thermal stability indicator
- **Antigenicity** - immunogenicity prediction
- **Flexibility/Rigidity** - structural domain boundaries
- **Transmembrane Regions** - membrane protein prediction
- **Disulfide Bonds** - structural stability
- **PTM Sites** - glycosylation, phosphorylation
- **Amino Acid Composition** - detailed breakdown

---

## PDB Structure Analysis / PDB结构分析

### Single File Analysis / 单文件分析

1. Click **"Select PDB File"**
2. Choose a `.pdb` file
3. Click **"Start Analysis"**

### Batch Analysis / 批量分析

1. Click **"Select Directory"**
2. Choose a folder containing PDB files
3. Click **"Start Analysis"**

### Analysis Results / 分析结果

- **Overall Score** - weighted geometric quality score
- **Chain Information** - for each protein chain:
  - Length
  - Geometric score
  - Clash detection (severe/mild)
  - Radius of gyration
  - Linker detection

### Ramachandran Plot / Ramachandran图

After PDB analysis, click **"Ramachandran Plot"** to generate.

---

## Comprehensive Analysis / 综合分析

Provides complete analysis combining both sequence and structure data:

1. Select PDB file and/or enter sequence
2. Click **"Generate Dashboard"** for summary
3. Click **"Generate All Charts"** for detailed plots

### Generated Charts / 生成的图表

- Hydrophobicity profile
- Flexibility profile
- Antigenicity profile
- PTM sites distribution
- Amino acid composition
- Comprehensive dashboard

---

## Results Export / 结果导出

### Export Formats / 导出格式

- **JSON** - Complete data for programmatic use
- **CSV** - Tabular data for Excel/spreadsheets
- **TXT** - Human-readable report

### Export Steps / 导出步骤

1. Go to **"Results View"** tab
2. Click **"Export Results"**
3. Choose format (JSON/CSV/TXT)
4. File saved to output folder

### Open Output Folder / 打开输出文件夹

Click **"Open Folder"** to view all generated files.

---

## Scoring System / 评分系统

### Score Ranges / 评分范围

| Score | Grade | Description |
|-------|-------|-------------|
| 80-100 | Excellent | High quality, ready for experimental validation |
| 60-79 | Good | Acceptable quality, minor issues possible |
| 40-59 | Average | Some concerns, careful review needed |
| 0-39 | Poor | Significant issues, redesign recommended |

### Score Factors / 评分因素

- Sequence length
- Isoelectric point
- Hydrophobicity
- Strong hydrophobic patches (penalty)
- Instability index (II>40 penalty)
- Charge balance
- PTM site potential (bonus)
- Disulfide bonds (bonus)

---

## Output Directory / 输出目录

Results are automatically saved to:

**Windows:**
```
C:\Users\[Username]\ProteinDesignEvaluator_Results\
```

**Linux/macOS:**
```
~/ProteinDesignEvaluator_Results/
```

### File Types / 文件类型

```
ProteinDesignEvaluator_Results/
├── figures/
│   ├── hydrophobicity_profile.png
│   ├── flexibility_profile.png
│   ├── antigenicity_profile.png
│   ├── ptm_sites.png
│   ├── aa_composition.png
│   ├── ramachandran.png
│   └── analysis_dashboard.png
├── analysis_20260611_143052.json
├── analysis_20260611_143052.csv
└── analysis_20260611_143052.txt
```

---

## FAQ / 常见问题

### Q: What amino acids are supported? / 支持哪些氨基酸?

A: 20 standard amino acids: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V

### Q: Minimum sequence length? / 最小序列长度?

A: At least 3 amino acids required for analysis.

### Q: PDB file format requirements? / PDB文件格式要求?

A: Standard PDB format with ATOM records. HETATM records supported.

### Q: How is the instability index calculated? / 不稳定指数如何计算?

A: Using Guruprasad method. II > 40 indicates potentially unstable protein.

### Q: What does antigenicity prediction mean? / 抗原性预测是什么意思?

A: Based on Kolaskar & Tongaonkar method. Higher values indicate stronger immunogenicity potential.

### Q: Can I use non-standard amino acids? / 可以使用非标准氨基酸吗?

A: Non-standard amino acids are automatically filtered out during analysis.

### Q: How are transmembrane regions predicted? / 如何预测跨膜区?

A: Using Kyte-Doolittle hydrophobicity with sliding window (window=20).

### Q: What is the aliphatic index? / 什么是脂肪族指数?

A: Measures thermal stability. Values > 100 indicate high thermal stability.

---

## Technical Details / 技术细节

### Algorithms Used / 使用的算法

| Feature | Algorithm/Method |
|---------|-----------------|
| pI | Iterative charge calculation |
| Hydrophobicity | Kyte-Doolittle scale |
| Antigenicity | Kolaskar & Tongaonkar |
| Flexibility | Karplus & Schulz |
| Instability | Guruprasad |
| Transmembrane | Kyte-Doolittle sliding window |
| N-Glycosylation | Asn-X-Ser/Thr motif |
| O-Glycosylation | Ser/Thr context analysis |
| Phosphorylation | S/T/Y with context scoring |

### References / 参考文献

1. Kyte J, Doolittle RF (1982). J Mol Biol.
2. Kolaskar AS, Tongaonkar PC (1990). FEBS Letters.
3. Guruprasad K et al. (1990). Protein Engineering.
4. Karplus PA, Schulz GE (1987). Naturwissenschaften.

---

**Document Version**: 3.0
**Last Updated**: 2026-06-11
