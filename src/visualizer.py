# -*- coding: utf-8 -*-
"""
Visualization Module v3.0
可视化模块

生成专业级的蛋白质分析图表
"""

import os
import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['axes.facecolor'] = '#f8f9fa'


def get_timestamp():
    """获取时间戳字符串"""
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


class ProteinVisualizer:
    """蛋白质分析可视化器 v3.0"""

    def __init__(self, output_dir):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.fig_dir = os.path.join(output_dir, "figures")
        os.makedirs(self.fig_dir, exist_ok=True)
        self.generated_files = []  # 记录生成的文件

    def plot_hydrophobicity_profile(self, sequence, profile, hydrophobic_patches=None,
                                    title="Hydrophobicity Profile", output_name=None):
        """绘制亲疏水性分布图"""
        length = len(profile)
        if length == 0:
            return None

        if output_name is None:
            output_name = f"hydrophobicity_profile_{get_timestamp()}.png"
        elif not output_name.endswith('.png'):
            output_name = f"{output_name}_{get_timestamp()}.png"
        else:
            name_part = output_name.replace('.png', '')
            output_name = f"{name_part}_{get_timestamp()}.png"

        fig, ax = plt.subplots(figsize=(max(12, length / 20), 5))

        positions = np.arange(length)
        ax.plot(positions, profile, color='#2196F3', linewidth=1.5, alpha=0.8)
        ax.fill_between(positions, profile, alpha=0.3, color='#2196F3')

        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax.axhline(y=1.6, color='red', linestyle='--', linewidth=1, alpha=0.5, label='TM threshold')
        ax.axhline(y=-1.6, color='green', linestyle='--', linewidth=1, alpha=0.5, label='Hydrophilic threshold')

        if hydrophobic_patches:
            for patch in hydrophobic_patches:
                start = patch['start']
                end = patch['end']
                rect = patches.Rectangle((start, min(profile)), end - start + 1,
                                        max(profile) - min(profile) + 2,
                                        linewidth=0, edgecolor='none',
                                        facecolor='red', alpha=0.2)
                ax.add_patch(rect)

        avg = np.mean(profile)
        ax.axhline(y=avg, color='orange', linestyle='-', linewidth=2, alpha=0.7, label=f'Mean: {avg:.2f}')

        ax.set_xlabel('Residue Position', fontsize=11)
        ax.set_ylabel('Hydrophobicity Index (Kyte-Doolittle)', fontsize=11)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.set_xlim(0, length)
        ax.set_ylim(min(profile) - 1, max(profile) + 2)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        # 添加时间戳避免覆盖
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_flexibility_profile(self, sequence, profile, flexible_regions=None,
                                  rigid_regions=None, title="Flexibility Profile",
                                  output_name="flexibility_profile.png"):
        """绘制柔性/刚性分布图"""
        length = len(profile)
        if length == 0:
            return None

        fig, ax = plt.subplots(figsize=(max(12, length / 20), 5))

        positions = np.arange(length)

        colors = []
        for val in profile:
            if val >= 0.70:
                colors.append('#F44336')
            elif val <= 0.50:
                colors.append('#4CAF50')
            else:
                colors.append('#FFC107')

        ax.bar(positions, profile, color=colors, alpha=0.8, width=0.8)

        ax.axhline(y=0.70, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Flexible threshold')
        ax.axhline(y=0.50, color='green', linestyle='--', linewidth=1, alpha=0.7, label='Rigid threshold')

        if flexible_regions:
            for region in flexible_regions:
                ax.axvspan(region['start'], region['end'], alpha=0.15, color='red')

        if rigid_regions:
            for region in rigid_regions:
                ax.axvspan(region['start'], region['end'], alpha=0.15, color='green')

        avg = np.mean(profile)
        ax.axhline(y=avg, color='blue', linestyle='-', linewidth=2, alpha=0.7, label=f'Mean: {avg:.3f}')

        ax.set_xlabel('Residue Position', fontsize=11)
        ax.set_ylabel('Flexibility Index', fontsize=11)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.set_xlim(0, length)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3, axis='y')

        legend_elements = [
            patches.Patch(facecolor='#F44336', label='Flexible'),
            patches.Patch(facecolor='#FFC107', label='Medium'),
            patches.Patch(facecolor='#4CAF50', label='Rigid'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

        plt.tight_layout()
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_antigenicity_profile(self, sequence, antigenicity_values, high_regions=None,
                                   title="Antigenicity Profile", output_name="antigenicity_profile.png"):
        """绘制抗原性分布图"""
        if len(antigenicity_values) == 0:
            return None

        window = 6
        length = len(sequence)
        positions = np.arange(window - 1, length)

        fig, ax = plt.subplots(figsize=(max(12, length / 20), 5))

        ax.plot(positions, antigenicity_values, color='#9C27B0', linewidth=1.5, alpha=0.8)
        ax.fill_between(positions, antigenicity_values, alpha=0.3, color='#9C27B0')

        ax.axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.7, label='High antigenicity (1.0)')
        ax.axhline(y=0.95, color='orange', linestyle='--', linewidth=1, alpha=0.7, label='Medium antigenicity (0.95)')

        if high_regions:
            for region in high_regions:
                ax.axvspan(region['start'], region['end'], alpha=0.2, color='red')

        avg = np.mean(antigenicity_values)
        ax.axhline(y=avg, color='blue', linestyle='-', linewidth=2, alpha=0.7, label=f'Mean: {avg:.3f}')

        ax.set_xlabel('Residue Position', fontsize=11)
        ax.set_ylabel('Antigenicity Index', fontsize=11)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.set_xlim(window - 1, length)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_property_radar(self, properties, title="Protein Property Radar",
                            output_name="property_radar.png"):
        """绘制综合理化性质雷达图"""
        categories = list(properties.keys())
        values = list(properties.values())

        values += values[:1]
        angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
        angles += angles[:1]

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

        ax.plot(angles, values, 'o-', linewidth=2, color='#2196F3')
        ax.fill(angles, values, alpha=0.25, color='#2196F3')

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=10)
        ax.set_ylim(0, 100)
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()
        output_path = os.path.join(self.fig_dir, output_name)
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_amino_acid_composition(self, composition, title="Amino Acid Composition",
                                    output_name="aa_composition.png"):
        """绘制氨基酸组成热图"""
        aa_order = list('ACDEFGHIKLMNPQRSTVWY')

        values = []
        labels = []
        for aa in aa_order:
            if aa in composition:
                values.append(composition[aa]['percentage'])
                labels.append(f"{aa}\n{composition[aa]['percentage']:.1f}%")
            else:
                values.append(0)
                labels.append(f"{aa}\n0%")

        fig, ax = plt.subplots(figsize=(14, 6))

        colors = []
        for v in values:
            if v > 9:
                colors.append('#F44336')
            elif v > 6:
                colors.append('#FF9800')
            elif v > 4:
                colors.append('#4CAF50')
            elif v > 2:
                colors.append('#FFC107')
            else:
                colors.append('#9E9E9E')

        bars = ax.bar(range(len(values)), values, color=colors, edgecolor='white', linewidth=1)

        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height,
                   f'{val:.1f}%', ha='center', va='bottom', fontsize=9)

        ax.set_xticks(range(len(aa_order)))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel('Percentage (%)', fontsize=11)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axhline(y=5, color='green', linestyle='--', alpha=0.5, label='Ideal mean (5%)')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_ramachandran(self, phi_psi_list, title="Ramachandran Plot",
                          output_name="ramachandran.png"):
        """绘制Ramachandran图"""
        fig, ax = plt.subplots(figsize=(10, 10))

        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xlabel('Phi (phi)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Psi (psi)', fontsize=14, fontweight='bold')
        ax.set_title(title, fontsize=16, fontweight='bold', pad=15)

        regions = [
            ([(-120, -120), (-120, -50), (-60, -50), (-60, -120)],
             '#32CD3233', '#228B22', 'Core alpha'),
            ([(-120, 50), (-120, 170), (-60, 170), (-60, 50)],
             '#4169E133', '#4169E1', 'Core beta'),
            ([(60, -60), (60, -120), (100, -120), (100, -60)],
             '#DAA52033', '#B8860B', 'L-alpha'),
            ([(-180, 50), (-100, 50), (-100, 180), (-180, 180)],
             '#9370DB33', '#9370DB', 'L-beta'),
        ]

        for coords, facecolor, edgecolor, label in regions:
            poly = Polygon(coords, facecolor=facecolor, edgecolor=edgecolor, linewidth=2)
            ax.add_patch(poly)
            cx = sum(c[0] for c in coords) / len(coords)
            cy = sum(c[1] for c in coords) / len(coords)
            ax.text(cx, cy, label, fontsize=11, ha='center', va='center', fontweight='bold')

        phi_vals = [p[0] for p in phi_psi_list if p[0] is not None and p[1] is not None]
        psi_vals = [p[1] for p in phi_psi_list if p[0] is not None and p[1] is not None]

        if phi_vals:
            positions = np.arange(len(phi_vals))
            scatter = ax.scatter(phi_vals, psi_vals, c=positions, cmap='plasma',
                              s=100, alpha=0.8, edgecolors='black', linewidths=1, zorder=10)

        ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.set_facecolor('#FAFAFA')

        plt.tight_layout()
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_ptm_sites(self, sequence, phosphorylation=None, n_glycosylation=None,
                       o_glycosylation=None, title="PTM Sites Distribution",
                       output_name="ptm_sites.png"):
        """绘制PTM位点分布图"""
        length = len(sequence)
        if length == 0:
            return None

        fig, ax = plt.subplots(figsize=(max(14, length / 15), 8))

        ax.set_xlim(0, length)
        ax.set_ylim(0, 5)

        for i, aa in enumerate(sequence):
            color = '#90A4AE'
            if aa in 'STY':
                if phosphorylation and any(s['position'] == i for s in phosphorylation.get('serine_threonine', []) + phosphorylation.get('tyrosine', [])):
                    color = '#F44336'
            if aa == 'N':
                if n_glycosylation and any(s['position'] == i for s in n_glycosylation.get('sites', [])):
                    color = '#2196F3'
            if aa in 'ST':
                if o_glycosylation and any(s['position'] == i for s in o_glycosylation.get('sites', [])):
                    color = '#4CAF50'

            ax.text(i + 0.5, 4.5, aa, ha='center', va='center', fontsize=6,
                   color=color, fontweight='bold')

        if phosphorylation:
            for site in phosphorylation.get('serine_threonine', []):
                ax.scatter(site['position'] + 0.5, 3, marker='s', s=100,
                          color='#F44336', edgecolors='black', linewidth=0.5, zorder=10)
            for site in phosphorylation.get('tyrosine', []):
                ax.scatter(site['position'] + 0.5, 3, marker='^', s=100,
                          color='#FF5722', edgecolors='black', linewidth=0.5, zorder=10)

        if n_glycosylation:
            for site in n_glycosylation.get('sites', []):
                ax.scatter(site['position'] + 0.5, 2, marker='o', s=150,
                          color='#2196F3', edgecolors='black', linewidth=1, zorder=10)
                ax.text(site['position'] + 0.5, 1.5, f"N-{site['position']}", ha='center', fontsize=7, color='#2196F3')

        if o_glycosylation:
            for site in o_glycosylation.get('sites', [])[:20]:
                ax.scatter(site['position'] + 0.5, 1, marker='D', s=80,
                          color='#4CAF50', edgecolors='black', linewidth=0.5, zorder=10)

        legend_elements = [
            patches.Patch(facecolor='#F44336', label=f"Phosphorylation (S/T: {len(phosphorylation.get('serine_threonine', []))}, Y: {len(phosphorylation.get('tyrosine', []))})"),
            patches.Patch(facecolor='#2196F3', label=f"N-Glycosylation ({len(n_glycosylation.get('sites', []))})"),
            patches.Patch(facecolor='#4CAF50', label=f"O-Glycosylation ({len(o_glycosylation.get('sites', []))})"),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

        ax.set_xticks([])
        ax.set_yticks([1, 2, 3, 4.5])
        ax.set_yticklabels(['O-Gly', 'N-Gly', 'Phospho', 'Sequence'], fontsize=9)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        plt.tight_layout()
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def create_analysis_dashboard(self, seq_results, pdb_results=None,
                                  output_name="analysis_dashboard.png"):
        """创建综合分析仪表板"""
        fig = plt.figure(figsize=(20, 16))

        fig.suptitle('Protein Comprehensive Analysis Report', fontsize=20, fontweight='bold', y=0.98)

        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)

        ax1 = fig.add_subplot(gs[0, 0])
        score = seq_results.get('comprehensive_score', {}).get('score', 0)
        grade = seq_results.get('comprehensive_score', {}).get('grade', 'N/A')
        colors = ['#4CAF50' if score >= 75 else '#FFC107' if score >= 50 else '#F44336']
        bars = ax1.barh(['Score'], [score], color=colors, height=0.5)
        ax1.set_xlim(0, 100)
        ax1.text(score + 2, 0, f'{score:.1f}', va='center', fontsize=16, fontweight='bold')
        ax1.set_title(f'Overall Score: {grade}', fontsize=12, fontweight='bold')
        ax1.axvline(x=75, color='green', linestyle='--', alpha=0.5)
        ax1.axvline(x=50, color='orange', linestyle='--', alpha=0.5)

        ax2 = fig.add_subplot(gs[0, 1])
        pI = seq_results.get('isoelectric_point', {}).get('pI', 7.0)
        ax2.barh(['pI'], [pI], color='#2196F3', height=0.5)
        ax2.set_xlim(0, 14)
        ax2.text(pI + 0.3, 0, f'pI = {pI}', va='center', fontsize=12, fontweight='bold')
        ax2.set_title('Isoelectric Point (pI)', fontsize=12, fontweight='bold')
        ax2.axvline(x=7.0, color='gray', linestyle='--', alpha=0.5)

        ax3 = fig.add_subplot(gs[0, 2])
        mw = seq_results.get('molecular_weight', {}).get('mw_kda', 0)
        ax3.barh(['MW'], [mw], color='#9C27B0', height=0.5)
        ax3.set_xlim(0, max(mw * 1.2, 100))
        ax3.text(mw + 0.5, 0, f'{mw:.1f} kDa', va='center', fontsize=12, fontweight='bold')
        ax3.set_title('Molecular Weight', fontsize=12, fontweight='bold')

        ax4 = fig.add_subplot(gs[0, 3])
        instability = seq_results.get('instability_index', {}).get('index', 0)
        is_stable = seq_results.get('instability_index', {}).get('is_stable', True)
        color = '#4CAF50' if is_stable else '#F44336'
        ax4.barh(['Stability'], [100 if is_stable else 0], color=color, height=0.5)
        ax4.set_xlim(0, 100)
        ax4.text(50 if is_stable else 50, 0, f'II={instability:.1f}', va='center',
                ha='center', fontsize=12, fontweight='bold')
        ax4.set_title('Instability Index', fontsize=12, fontweight='bold')

        ax5 = fig.add_subplot(gs[1, :2])
        profile = seq_results.get('hydrophobicity_profile', {}).get('profile', [])
        if profile:
            positions = np.arange(len(profile))
            ax5.plot(positions, profile, color='#2196F3', linewidth=1)
            ax5.fill_between(positions, profile, alpha=0.3, color='#2196F3')
            ax5.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            avg = np.mean(profile)
            ax5.axhline(y=avg, color='orange', linestyle='-', linewidth=2, label=f'Mean: {avg:.2f}')
            ax5.set_title('Hydrophobicity Profile', fontsize=12, fontweight='bold')
            ax5.set_xlabel('Residue Position')
            ax5.set_ylabel('Kyte-Doolittle Hydrophobicity')
            ax5.legend()

        ax6 = fig.add_subplot(gs[1, 2:])
        comp = seq_results.get('composition', {})
        aa_order = list('ACDEFGHIKLMNPQRSTVWY')
        values = [comp.get(aa, {}).get('percentage', 0) for aa in aa_order]
        colors = ['#4FC3F7' if v > 5 else '#90A4AE' for v in values]
        ax6.bar(aa_order, values, color=colors, edgecolor='white')
        ax6.axhline(y=5, color='red', linestyle='--', alpha=0.5, label='Ideal 5%')
        ax6.set_title('Amino Acid Composition', fontsize=12, fontweight='bold')
        ax6.set_ylabel('Percentage (%)')
        ax6.legend()

        ax7 = fig.add_subplot(gs[2, :2])
        n_glyco = len(seq_results.get('n_glycosylation', {}).get('sites', []))
        o_glyco = len(seq_results.get('o_glycosylation', {}).get('sites', []))
        phos_st = len(seq_results.get('phosphorylation', {}).get('serine_threonine', []))
        phos_y = len(seq_results.get('phosphorylation', {}).get('tyrosine', []))
        disul = seq_results.get('disulfide_bonds', {}).get('estimated_bonds', 0)

        categories = ['N-Gly', 'O-Gly', 'Phospho\n(S/T)', 'Phospho\n(Y)', 'Disulfide']
        values = [n_glyco, o_glyco, phos_st, phos_y, disul]
        colors = ['#2196F3', '#4CAF50', '#F44336', '#FF5722', '#FFC107']
        ax7.bar(categories, values, color=colors, edgecolor='white')
        for i, v in enumerate(values):
            ax7.text(i, v + 0.1, str(v), ha='center', fontsize=10, fontweight='bold')
        ax7.set_title('PTM Sites', fontsize=12, fontweight='bold')
        ax7.set_ylabel('Count')

        ax8 = fig.add_subplot(gs[2, 2:])
        metrics = [
            ('Hydrophobicity', seq_results.get('hydrophobicity', {}).get('average', 0), 2.0),
            ('Aliphatic Index', seq_results.get('aliphatic_index', {}).get('index', 0), 100),
            ('Antigenicity', seq_results.get('antigenicity', {}).get('score', 0), 1.0),
            ('Flexibility', seq_results.get('flexibility', {}).get('average', 0), 0.7),
        ]
        names = [m[0] for m in metrics]
        vals = [m[1] for m in metrics]
        maxes = [m[2] for m in metrics]
        normalized = [v / m * 100 for v, m in zip(vals, maxes)]
        ax8.barh(names, normalized, color=['#2196F3', '#9C27B0', '#FF9800', '#4CAF50'])
        ax8.set_xlim(0, 100)
        ax8.set_title('Physicochemical Properties (Normalized)', fontsize=12, fontweight='bold')
        for i, (n, v) in enumerate(zip(names, vals)):
            ax8.text(normalized[i] + 2, i, f'{v:.2f}', va='center', fontsize=9)

        ax9 = fig.add_subplot(gs[3, :])
        ax9.axis('off')

        findings = []
        if seq_results.get('transmembrane', {}).get('has_transmembrane'):
            regions = seq_results['transmembrane']['regions']
            findings.append(f"Predicted {len(regions)} transmembrane region(s)")
        if seq_results.get('disulfide_bonds', {}).get('estimated_bonds', 0) > 0:
            findings.append(f"Possible {seq_results['disulfide_bonds']['estimated_bonds']} disulfide bond(s)")
        if seq_results.get('n_glycosylation', {}).get('count', 0) > 0:
            findings.append(f"{seq_results['n_glycosylation']['count']} N-glycosylation site(s)")
        if seq_results.get('instability_index', {}).get('is_stable', True):
            findings.append("Protein predicted as STABLE")
        else:
            findings.append("Protein predicted as UNSTABLE")

        hydrophobic_patches = seq_results.get('hydrophobicity_profile', {}).get('hydrophobic_patches', [])
        if hydrophobic_patches:
            findings.append(f"Found {len(hydrophobic_patches)} strong hydrophobic patch(es)")

        if pdb_results:
            geo_score = pdb_results.get('overall_score', 0)
            findings.append(f"Structure score: {geo_score:.1f}/100")

        findings_text = " | ".join(findings) if findings else "No significant features found"
        ax9.text(0.5, 0.5, findings_text, ha='center', va='center', fontsize=12,
                bbox=dict(boxstyle='round', facecolor='#E3F2FD', edgecolor='#2196F3', linewidth=2))

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path

    def plot_batch_comparison(self, results, metric='comprehensive_score',
                              title="Batch Analysis Comparison", output_name="batch_comparison.png"):
        """绘制批量分析结果对比图"""
        if not results:
            return None

        n = min(len(results), 30)
        fig, axes = plt.subplots(2, 1, figsize=(14, 10))

        names = []
        scores = []
        for r in results[:n]:
            if isinstance(r, dict):
                name = r.get('id', r.get('basename', f'Sample {results.index(r)+1}'))[:30]
                if metric == 'comprehensive_score':
                    score = r.get('comprehensive_score', {}).get('score', 0)
                else:
                    score = r.get(metric, 0)
            else:
                name = f"Sample {results.index(r)+1}"
                score = 0
            names.append(name)
            scores.append(score)

        sorted_data = sorted(zip(names, scores), key=lambda x: x[1], reverse=True)
        names, scores = zip(*sorted_data)

        colors = ['#4FC3F7' if s >= 75 else '#81C784' if s >= 50 else '#FFC107' if s >= 30 else '#F44336'
                 for s in scores]

        axes[0].barh(range(len(names)), scores, color=colors, edgecolor='white')
        axes[0].set_yticks(range(len(names)))
        axes[0].set_yticklabels(names, fontsize=8)
        axes[0].set_xlabel('Score', fontsize=11)
        axes[0].set_title(title, fontsize=14, fontweight='bold')
        axes[0].invert_yaxis()
        axes[0].set_xlim(0, 100)

        for i, s in enumerate(scores):
            axes[0].text(s + 1, i, f'{s:.1f}', va='center', fontsize=8)

        axes[1].hist(scores, bins=15, color='#2196F3', edgecolor='white', alpha=0.8)
        axes[1].axvline(np.mean(scores), color='red', linestyle='--', linewidth=2,
                        label=f'Mean: {np.mean(scores):.1f}')
        axes[1].axvline(np.median(scores), color='orange', linestyle='--', linewidth=2,
                        label=f'Median: {np.median(scores):.1f}')
        axes[1].set_xlabel('Score', fontsize=11)
        axes[1].set_ylabel('Count', fontsize=11)
        axes[1].set_title('Score Distribution', fontsize=14, fontweight='bold')
        axes[1].legend()

        plt.tight_layout()
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        name_base = output_name.replace('.png', '')
        output_path = os.path.join(self.fig_dir, f"{name_base}_{timestamp}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return output_path
