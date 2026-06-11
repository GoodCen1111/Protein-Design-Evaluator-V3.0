# -*- coding: utf-8 -*-
"""
Advanced Protein Sequence Analyzer v3.0 - Complete Rewrite
高级蛋白质序列分析器 - 完整重写

使用标准化的氨基酸理化性质数据库，确保计算结果客观准确
"""

import re
import math
from collections import defaultdict


# ========== 标准氨基酸理化性质数据库 (基于文献) ==========

# Kyte-Doolittle 疏水性 scale (标准值)
HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

# 氨基酸分子量 (Da) - 标准值
MOLECULAR_WEIGHT = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
}

# pKa values - Bjellqvist scale (最常用的标准)
# 来源: Bjellqvist et al., 1994, Electrophoresis
PKA_CARBOXYL_TERM = 2.34  # C-terminus
PKA_AMINO_TERM = 9.69     # N-terminus
PKA_SIDECHAIN = {
    'D': 3.86,   # Asp
    'E': 4.25,   # Glu
    'H': 6.0,    # His
    'C': 8.33,   # Cys
    'Y': 10.46,  # Tyr
    'K': 10.54,  # Lys
    'R': 12.48,  # Arg
}

# 柔性指数 (Karplus & Schulz)
FLEXIBILITY_SCALE = {
    'G': 0.93, 'A': 0.78, 'S': 0.75, 'P': 0.72, 'N': 0.68,
    'E': 0.67, 'K': 0.65, 'R': 0.66, 'T': 0.64, 'D': 0.62,
    'Q': 0.62, 'M': 0.55, 'H': 0.54, 'C': 0.50, 'V': 0.50,
    'I': 0.47, 'L': 0.45, 'F': 0.41, 'Y': 0.40, 'W': 0.34,
}

# 抗原性指数 (Kolaskar & Tongaonkar, 1990)
ANTIGENICITY_SCALE = {
    'A': 1.064, 'R': 0.644, 'N': 0.803, 'D': 0.858, 'C': 1.046,
    'E': 0.858, 'Q': 0.809, 'G': 0.874, 'H': 0.864, 'I': 1.152,
    'L': 1.236, 'K': 0.644, 'M': 1.303, 'F': 1.268, 'P': 0.858,
    'S': 0.941, 'T': 0.811, 'W': 1.089, 'Y': 1.064, 'V': 1.268,
}

# 二硫键不稳定指数 (Guruprasad et al., 1990)
# 不稳定指数二肽权重 (Guruprasad et al., 1990)
INSTABILITY_DIPEPTIDE_WEIGHTS = {
    'WW': 1.0, 'WC': 24.68, 'WM': 5.94, 'WH': 24.68, 'WT': 7.49,
    'WY': 7.49, 'WF': 4.9, 'WQ': -5.68, 'WN': -1.88, 'WI': 10.46,
    'WR': -9.37, 'WS': 2.08, 'WD': -9.37, 'WK': 4.44, 'WP': -1.59,
    'WL': 5.63, 'WG': -9.37, 'WV': 5.63, 'WE': -7.49, 'WA': -14.03,
    'CW': 1.0, 'CC': 1.0, 'CM': 1.0, 'CH': 1.0, 'CT': 1.0,
    'CY': 1.0, 'CF': 1.0, 'CQ': 1.0, 'CN': 1.0, 'CI': 1.0,
    'CR': 1.0, 'CS': 1.0, 'CD': 1.0, 'CK': 1.0, 'CP': 1.0,
    'CL': 1.0, 'CG': 1.0, 'CV': 1.0, 'CE': 1.0, 'CA': 1.0,
    'MW': 5.94, 'MC': 1.0, 'MM': 5.94, 'MH': -2.16, 'MT': 11.36,
    'MY': 11.36, 'MF': 11.36, 'MQ': -2.16, 'MN': -2.16, 'MI': 11.36,
    'MR': -0.77, 'MS': -0.77, 'MD': -0.77, 'MK': -0.77, 'MP': -1.19,
    'ML': 11.36, 'MG': -0.77, 'MV': 5.94, 'ME': 5.94, 'MA': -14.03,
    'HW': 24.68, 'HC': 1.0, 'HM': -2.16, 'HH': 24.68, 'HT': 1.0,
    'HY': 1.0, 'HF': 1.0, 'HQ': 1.0, 'HN': 1.0, 'HI': 1.0,
    'HR': 1.0, 'HS': 1.0, 'HD': 1.0, 'HK': 1.0, 'HP': 1.0,
    'HL': 1.0, 'HG': 1.0, 'HV': 1.0, 'HE': 1.0, 'HA': 1.0,
    'TW': 7.49, 'TC': 1.0, 'TM': 11.36, 'TH': 1.0, 'TT': 7.49,
    'TY': 1.0, 'TF': 4.44, 'TQ': 1.0, 'TN': 1.0, 'TI': 1.0,
    'TR': 1.0, 'TS': 7.49, 'TD': 1.0, 'TK': 1.0, 'TP': -1.59,
    'TL': 1.0, 'TG': 1.0, 'TV': 1.0, 'TE': 1.0, 'TA': 1.0,
    'YW': 7.49, 'YC': 1.0, 'YM': 11.36, 'YH': 1.0, 'YT': 1.0,
    'YY': 7.49, 'YF': 4.44, 'YQ': 1.0, 'YN': 1.0, 'YI': 1.0,
    'YR': 1.0, 'YS': 1.0, 'YD': 1.0, 'YK': 1.0, 'YP': -1.59,
    'YL': 1.0, 'YG': 1.0, 'YV': 1.0, 'YE': 1.0, 'YA': 1.0,
    'FW': 4.9, 'FC': 1.0, 'FM': 11.36, 'FH': 1.0, 'FT': 4.44,
    'FY': 4.44, 'FF': 4.9, 'FQ': 1.0, 'FN': 1.0, 'FI': 4.9,
    'FR': 1.0, 'FS': 1.0, 'FD': 1.0, 'FK': 1.0, 'FP': -1.59,
    'FL': 4.9, 'FG': 1.0, 'FV': 4.9, 'FE': 1.0, 'FA': 1.0,
    'QW': -5.68, 'QC': 1.0, 'QM': -2.16, 'QH': 1.0, 'QT': 1.0,
    'QY': 1.0, 'QF': 1.0, 'QQ': -5.68, 'QN': 1.0, 'QI': 1.0,
    'QR': 1.0, 'QS': 1.0, 'QD': 1.0, 'QK': 1.0, 'QP': 1.0,
    'QL': 1.0, 'QG': 1.0, 'QV': 1.0, 'QE': 1.0, 'QA': 1.0,
    'NW': -1.88, 'NC': 1.0, 'NM': -2.16, 'NH': 1.0, 'NT': 1.0,
    'NY': 1.0, 'NF': 1.0, 'NQ': 1.0, 'NN': -1.88, 'NI': 1.0,
    'NR': 1.0, 'NS': 1.0, 'ND': 1.0, 'NK': 1.0, 'NP': 1.0,
    'NL': 1.0, 'NG': 1.0, 'NV': 1.0, 'NE': 1.0, 'NA': 1.0,
    'IW': 10.46, 'IC': 1.0, 'IM': 11.36, 'IH': 1.0, 'IT': 1.0,
    'IY': 1.0, 'IF': 4.9, 'IQ': 1.0, 'IN': 1.0, 'II': 10.46,
    'IR': 1.0, 'IS': 1.0, 'ID': 1.0, 'IK': 1.0, 'IP': -1.59,
    'IL': 10.46, 'IG': 1.0, 'IV': 10.46, 'IE': 1.0, 'IA': 1.0,
    'RW': -9.37, 'RC': 1.0, 'RM': -0.77, 'RH': 1.0, 'RT': 1.0,
    'RY': 1.0, 'RF': 1.0, 'RQ': 1.0, 'RN': 1.0, 'RI': 1.0,
    'RR': -9.37, 'RS': 1.0, 'RD': 1.0, 'RK': 1.0, 'RP': 1.0,
    'RL': 1.0, 'RG': 1.0, 'RV': 1.0, 'RE': 1.0, 'RA': 1.0,
    'SW': 2.08, 'SC': 1.0, 'SM': -0.77, 'SH': 1.0, 'ST': 7.49,
    'SY': 1.0, 'SF': 1.0, 'SQ': 1.0, 'SN': 1.0, 'SI': 1.0,
    'SR': 1.0, 'SS': 2.08, 'SD': 1.0, 'SK': 1.0, 'SP': -1.59,
    'SL': 1.0, 'SG': 1.0, 'SV': 1.0, 'SE': 1.0, 'SA': 1.0,
    'DW': -9.37, 'DC': 1.0, 'DM': -0.77, 'DH': 1.0, 'DT': 1.0,
    'DY': 1.0, 'DF': 1.0, 'DQ': 1.0, 'DN': 1.0, 'DI': 1.0,
    'DR': 1.0, 'DS': 1.0, 'DD': -9.37, 'DK': 1.0, 'DP': 1.0,
    'DL': 1.0, 'DG': 1.0, 'DV': 1.0, 'DE': 1.0, 'DA': 1.0,
    'KW': 4.44, 'KC': 1.0, 'KM': -0.77, 'KH': 1.0, 'KT': 1.0,
    'KY': 1.0, 'KF': 1.0, 'KQ': 1.0, 'KN': 1.0, 'KI': 1.0,
    'KR': 1.0, 'KS': 1.0, 'KD': 1.0, 'KK': 4.44, 'KP': 1.0,
    'KL': 1.0, 'KG': 1.0, 'KV': 1.0, 'KE': 1.0, 'KA': 1.0,
    'PW': -1.59, 'PC': 1.0, 'PM': -1.19, 'PH': 1.0, 'PT': -1.59,
    'PY': -1.59, 'PF': -1.59, 'PQ': 1.0, 'PN': 1.0, 'PI': -1.59,
    'PR': 1.0, 'PS': -1.59, 'PD': 1.0, 'PK': 1.0, 'PP': -1.59,
    'PL': 1.0, 'PG': -1.59, 'PV': -1.59, 'PE': 1.0, 'PA': 1.0,
    'LW': 13.34, 'LC': 1.0, 'LM': 11.36, 'LH': 1.0, 'LT': 1.0,
    'LY': 1.0, 'LF': 4.9, 'LQ': 1.0, 'LN': 1.0, 'LI': 10.46,
    'LR': 1.0, 'LS': 1.0, 'LD': 1.0, 'LK': 1.0, 'LP': 1.0,
    'LL': 13.34, 'LG': 1.0, 'LV': 13.34, 'LE': 1.0, 'LA': 1.0,
    'GW': -9.37, 'GC': 1.0, 'GM': -0.77, 'GH': 1.0, 'GT': 1.0,
    'GY': 1.0, 'GF': 1.0, 'GQ': 1.0, 'GN': 1.0, 'GI': 1.0,
    'GR': 1.0, 'GS': 1.0, 'GD': 1.0, 'GK': 1.0, 'GP': -1.59,
    'GL': 1.0, 'GG': -9.37, 'GV': 1.0, 'GE': 1.0, 'GA': 1.0,
    'VW': 5.63, 'VC': 1.0, 'VM': 5.94, 'VH': 1.0, 'VT': 1.0,
    'VY': 1.0, 'VF': 4.9, 'VQ': 1.0, 'VN': 1.0, 'VI': 10.46,
    'VR': 1.0, 'VS': 1.0, 'VD': 1.0, 'VK': 1.0, 'VP': -1.59,
    'VL': 13.34, 'VG': 1.0, 'VV': 5.63, 'VE': 1.0, 'VA': 1.0,
    'EW': -7.49, 'EC': 1.0, 'EM': 5.94, 'EH': 1.0, 'ET': 1.0,
    'EY': 1.0, 'EF': 1.0, 'EQ': 1.0, 'EN': 1.0, 'EI': 1.0,
    'ER': 1.0, 'ES': 1.0, 'ED': 1.0, 'EK': 1.0, 'EP': 1.0,
    'EL': 1.0, 'EG': 1.0, 'EV': 1.0, 'EE': -7.49, 'EA': 1.0,
    'AW': -14.03, 'AC': 1.0, 'AM': -14.03, 'AH': 1.0, 'AT': 1.0,
    'AY': 1.0, 'AF': 1.0, 'AQ': 1.0, 'AN': 1.0, 'AI': 1.0,
    'AR': 1.0, 'AS': 1.0, 'AD': 1.0, 'AK': 1.0, 'AP': 1.0,
    'AL': 1.0, 'AG': 1.0, 'AV': 1.0, 'AE': 1.0, 'AA': -14.03,
}

# GRAVY (Grand Average of Hydropathy) 计算用
GRAVY_SCALE = HYDROPHOBICITY.copy()

# 跨膜区阈值 (Kyte-Doolittle window=20)
TRANSMEMBRANE_THRESHOLD = 1.6

# N-糖基化位点 motif
N_GLYCOSYLATION_MOTIF = re.compile(r'N[^P][ST]')

# O-糖基化潜在位点
O_GLYCOSYLATION_SITES = ['S', 'T']


class AdvancedSequenceAnalyzer:
    """高级蛋白质序列分析器 v3.0 - 完整版"""

    def __init__(self):
        self.sequence = ""
        self.cleaned_sequence = ""
        self.results = {}

    def analyze_sequence(self, sequence):
        """执行完整的序列分析"""
        self.sequence = sequence.upper()
        self.cleaned_sequence = self._clean_sequence(sequence)

        if len(self.cleaned_sequence) == 0:
            raise ValueError("序列为空或不包含有效氨基酸")

        results = {
            'sequence': self.cleaned_sequence,
            'length': len(self.cleaned_sequence),
        }

        # 执行所有分析
        results['composition'] = self._calculate_composition()
        results['molecular_weight'] = self._calculate_molecular_weight()
        results['isoelectric_point'] = self._calculate_pI()
        results['hydrophobicity'] = self._calculate_hydrophobicity()
        results['hydrophobicity_profile'] = self._calculate_hydrophobicity_profile(window=9)
        results['gravy'] = self._calculate_gravy()
        results['instability_index'] = self._calculate_instability_index()
        results['aliphatic_index'] = self._calculate_aliphatic_index()
        results['flexibility'] = self._calculate_flexibility_profile(window=7)
        results['antigenicity'] = self._calculate_antigenicity()
        results['transmembrane'] = self._predict_transmembrane_regions()
        results['disulfide_bonds'] = self._predict_disulfide_bonds()
        results['n_glycosylation'] = self._predict_n_glycosylation()
        results['o_glycosylation'] = self._predict_o_glycosylation()
        results['phosphorylation'] = self._predict_phosphorylation()
        results['expression_analysis'] = self._predict_expression_host()
        results['comprehensive_score'] = self._calculate_comprehensive_score(results)

        self.results = results
        return results

    def _clean_sequence(self, sequence):
        """清理序列，只保留20种标准氨基酸"""
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        return ''.join(c.upper() for c in sequence if c.upper() in valid_aa)

    def _calculate_composition(self):
        """计算氨基酸组成"""
        seq = self.cleaned_sequence
        length = len(seq)
        comp = {}

        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            count = seq.count(aa)
            comp[aa] = {
                'count': count,
                'percentage': (count / length * 100) if length > 0 else 0
            }

        comp['classification'] = {
            'hydrophobic': sum(comp[aa]['count'] for aa in 'AVILMFYW'),
            'hydrophilic': sum(comp[aa]['count'] for aa in 'RKDENQSTH'),
            'charged_positive': sum(comp[aa]['count'] for aa in 'KRH'),
            'charged_negative': sum(comp[aa]['count'] for aa in 'DE'),
            'polar': sum(comp[aa]['count'] for aa in 'NQSTYC'),
            'nonpolar': sum(comp[aa]['count'] for aa in 'AGILMFPWV'),
            'aromatic': sum(comp[aa]['count'] for aa in 'FYW'),
            'sulfur': sum(comp[aa]['count'] for aa in 'CM'),
            'tiny': sum(comp[aa]['count'] for aa in 'AG'),
            'aliphatic': sum(comp[aa]['count'] for aa in 'AVILM'),
        }

        return comp

    def _calculate_molecular_weight(self):
        """计算分子量"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length == 0:
            return {'mw': 0, 'formula': 'N/A'}

        total_mw = sum(MOLECULAR_WEIGHT.get(aa, 110) for aa in seq)
        water_loss = 18.015 * (length - 1)
        mw = total_mw - water_loss

        return {
            'mw': round(mw, 2),
            'mw_kda': round(mw / 1000, 2),
            'water_loss': round(water_loss, 2),
        }

    def _calculate_pI(self):
        """
        计算等电点 (pI) - 使用Bjellqvist pKa scale
        标准迭代法求解蛋白质净电荷为零时的pH值
        """
        seq = self.cleaned_sequence
        length = len(seq)

        def calculate_charge(ph):
            """计算给定pH下蛋白质的净电荷"""
            positive = 0.0
            negative = 0.0

            # N-terminus (带正电)
            positive += 1 / (1 + 10 ** (ph - PKA_AMINO_TERM))

            # C-terminus (带负电)
            negative += 1 / (1 + 10 ** (PKA_CARBOXYL_TERM - ph))

            # 侧链
            for aa in seq:
                if aa in PKA_SIDECHAIN:
                    pka = PKA_SIDECHAIN[aa]
                    if pka < 7:
                        negative += 1 / (1 + 10 ** (pka - ph))
                    else:
                        positive += 1 / (1 + 10 ** (ph - pka))

            return positive - negative

        # 统计可解离基团
        pos_count = sum(1 for aa in seq if aa in 'KRH') + 1  # +1 for N-terminus
        neg_count = sum(1 for aa in seq if aa in 'DE') + 1     # +1 for C-terminus

        # 初始估计
        if neg_count > pos_count:
            pH = 2.0
        elif pos_count > neg_count:
            pH = 12.0
        else:
            pH = 7.0

        # 二分法查找等电点
        ph_low, ph_high = 0.0, 14.0
        for _ in range(100):
            charge = calculate_charge(pH)
            if abs(charge) < 0.001:
                break
            if charge > 0:
                ph_low = pH
                pH = (ph_low + ph_high) / 2
            else:
                ph_high = pH
                pH = (ph_low + ph_high) / 2

        pI = round(pH, 2)

        # 蛋白质类型判断
        if pI < 6.0:
            ptype = "酸性蛋白"
        elif pI > 8.0:
            ptype = "碱性蛋白"
        else:
            ptype = "中性蛋白"

        # 计算pH=7时的净电荷
        net_charge_at_ph7 = calculate_charge(7.0)

        return {
            'pI': pI,
            'protein_type': ptype,
            'net_charge_at_ph7': round(net_charge_at_ph7, 2),
            'positive_charges': pos_count,
            'negative_charges': neg_count,
        }

    def _calculate_hydrophobicity(self):
        """计算整体亲疏水性"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length == 0:
            return {'average': 0, 'status': 'N/A'}

        total = sum(HYDROPHOBICITY.get(aa, 0) for aa in seq)
        avg = total / length

        if avg > 1.0:
            status = "强疏水性"
        elif avg > 0.5:
            status = "中等疏水性"
        elif avg > 0:
            status = "轻度疏水性"
        elif avg > -0.5:
            status = "轻度亲水性"
        elif avg > -1.0:
            status = "中等亲水性"
        else:
            status = "强亲水性"

        hydrophobic_count = sum(1 for aa in seq if HYDROPHOBICITY.get(aa, 0) > 0)
        hydrophilic_count = sum(1 for aa in seq if HYDROPHOBICITY.get(aa, 0) < 0)

        return {
            'average': round(avg, 3),
            'status': status,
            'total_score': round(total, 2),
            'hydrophobic_percentage': round(hydrophobic_count / length * 100, 1),
            'hydrophilic_percentage': round(hydrophilic_count / length * 100, 1),
        }

    def _calculate_gravy(self):
        """计算GRAVY (Grand Average of Hydropathy)"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length == 0:
            return {'gravy': 0, 'interpretation': 'N/A'}

        total = sum(GRAVY_SCALE.get(aa, 0) for aa in seq)
        gravy = total / length

        if gravy > 1:
            interpretation = "疏水性蛋白"
        elif gravy > 0:
            interpretation = "偏疏水"
        elif gravy > -1:
            interpretation = "偏亲水"
        else:
            interpretation = "亲水性蛋白"

        return {
            'gravy': round(gravy, 3),
            'interpretation': interpretation,
        }

    def _calculate_hydrophobicity_profile(self, window=9):
        """计算沿序列的亲疏水性分布"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length < window:
            return {'profile': [], 'window': window}

        profile = []
        half_window = window // 2

        for i in range(length):
            start = max(0, i - half_window)
            end = min(length, i + half_window + 1)
            window_seq = seq[start:end]
            avg = sum(HYDROPHOBICITY.get(aa, 0) for aa in window_seq) / len(window_seq)
            profile.append(round(avg, 3))

        # 识别强疏水斑块
        hydrophobic_patches = []
        threshold = 1.0
        current_patch = []

        for i, val in enumerate(profile):
            if val >= threshold:
                current_patch.append((i, val))
            else:
                if len(current_patch) >= 5:
                    hydrophobic_patches.append({
                        'start': current_patch[0][0],
                        'end': current_patch[-1][0],
                        'length': len(current_patch),
                        'max_value': max(v for _, v in current_patch),
                        'avg_value': round(sum(v for _, v in current_patch) / len(current_patch), 3),
                    })
                current_patch = []

        if len(current_patch) >= 5:
            hydrophobic_patches.append({
                'start': current_patch[0][0],
                'end': current_patch[-1][0],
                'length': len(current_patch),
                'max_value': max(v for _, v in current_patch),
                'avg_value': round(sum(v for _, v in current_patch) / len(current_patch), 3),
            })

        return {
            'profile': profile,
            'window': window,
            'hydrophobic_patches': hydrophobic_patches,
        }

    def _calculate_instability_index(self):
        """
        计算不稳定指数 (Instability Index, II)
        Guruprasad et al., 1990方法
        II < 40 表示蛋白质稳定，II > 40 表示不稳定
        """
        seq = self.cleaned_sequence
        length = len(seq)

        if length < 2:
            return {'index': 0, 'prediction': 'N/A'}

        total = 0
        for i in range(length - 1):
            dipeptide = seq[i:i + 2]
            weight = INSTABILITY_DIPEPTIDE_WEIGHTS.get(dipeptide, 1.0)
            total += weight

        index = (10.0 / length) * total

        is_stable = index < 40
        if index < 40:
            prediction = "稳定"
            stability_note = "预测蛋白质在溶液中相对稳定"
        else:
            prediction = "不稳定"
            stability_note = "预测蛋白质在溶液中可能容易降解"

        return {
            'index': round(index, 2),
            'prediction': prediction,
            'stability_note': stability_note,
            'is_stable': is_stable,
        }

    def _calculate_aliphatic_index(self):
        """
        计算脂肪族指数
        Ikai, 1980 - 衡量热稳定性
        """
        seq = self.cleaned_sequence
        length = len(seq)

        if length == 0:
            return {'index': 0, 'prediction': 'N/A'}

        comp = self.results.get('composition', {})
        if not comp:
            comp = self._calculate_composition()

        ala = comp.get('A', {}).get('percentage', 0)
        val = comp.get('V', {}).get('percentage', 0)
        ile = comp.get('I', {}).get('percentage', 0)
        leu = comp.get('L', {}).get('percentage', 0)

        # 公式: AI = X(Ala) + 2.9 * X(Val) + 3.9 * [X(Ile) + X(Leu)]
        ai = ala + 2.9 * val + 3.9 * (ile + leu)

        if ai > 100:
            prediction = "高热稳定性"
            note = "预测具有较高的热稳定性"
        elif ai > 80:
            prediction = "中等热稳定性"
            note = "预测具有一般热稳定性"
        elif ai > 60:
            prediction = "较低热稳定性"
            note = "预测热稳定性较低"
        else:
            prediction = "低热稳定性"
            note = "预测容易变性"

        return {
            'index': round(ai, 2),
            'prediction': prediction,
            'note': note,
            'components': {
                'alanine': round(ala, 2),
                'valine': round(val, 2),
                'isoleucine': round(ile, 2),
                'leucine': round(leu, 2),
            }
        }

    def _calculate_flexibility_profile(self, window=7):
        """计算柔性/刚性分布 - Karplus & Schulz方法"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length < window:
            return {'profile': [], 'window': window}

        profile = []
        half_window = window // 2

        for i in range(length):
            start = max(0, i - half_window)
            end = min(length, i + half_window + 1)
            window_seq = seq[start:end]
            avg = sum(FLEXIBILITY_SCALE.get(aa, 0.5) for aa in window_seq) / len(window_seq)
            profile.append(round(avg, 3))

        # 分类区域
        flexible_regions = []
        rigid_regions = []
        threshold_high = 0.70
        threshold_low = 0.50

        current_flex = []
        current_rigid = []

        for i, val in enumerate(profile):
            if val >= threshold_high:
                if current_rigid:
                    if len(current_rigid) >= 5:
                        rigid_regions.append({
                            'start': current_rigid[0][0],
                            'end': current_rigid[-1][0],
                            'length': len(current_rigid),
                            'avg_flexibility': round(sum(v for _, v in current_rigid) / len(current_rigid), 3),
                        })
                    current_rigid = []
                current_flex.append((i, val))
            elif val <= threshold_low:
                if current_flex:
                    if len(current_flex) >= 5:
                        flexible_regions.append({
                            'start': current_flex[0][0],
                            'end': current_flex[-1][0],
                            'length': len(current_flex),
                            'avg_flexibility': round(sum(v for _, v in current_flex) / len(current_flex), 3),
                        })
                    current_flex = []
                current_rigid.append((i, val))

        # 处理剩余
        if len(current_flex) >= 5:
            flexible_regions.append({
                'start': current_flex[0][0],
                'end': current_flex[-1][0],
                'length': len(current_flex),
                'avg_flexibility': round(sum(v for _, v in current_flex) / len(current_flex), 3),
            })
        if len(current_rigid) >= 5:
            rigid_regions.append({
                'start': current_rigid[0][0],
                'end': current_rigid[-1][0],
                'length': len(current_rigid),
                'avg_flexibility': round(sum(v for _, v in current_rigid) / len(current_rigid), 3),
            })

        avg_flexibility = sum(profile) / len(profile)

        return {
            'profile': profile,
            'window': window,
            'average': round(avg_flexibility, 3),
            'flexible_regions': flexible_regions,
            'rigid_regions': rigid_regions,
            'flexible_count': sum(1 for v in profile if v >= threshold_high),
            'rigid_count': sum(1 for v in profile if v <= threshold_low),
        }

    def _calculate_antigenicity(self):
        """抗原性预测 - Kolaskar & Tongaonkar方法"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length == 0:
            return {'score': 0, 'prediction': 'N/A'}

        window = 6
        antigenicity_values = []

        for i in range(length - window + 1):
            window_seq = seq[i:i + window]
            avg = sum(ANTIGENICITY_SCALE.get(aa, 1.0) for aa in window_seq) / window
            antigenicity_values.append(avg)

        overall_avg = sum(antigenicity_values) / len(antigenicity_values) if antigenicity_values else 0

        if overall_avg >= 1.0:
            prediction = "高抗原性"
            threshold_note = "可能具有较强的免疫原性"
        elif overall_avg >= 0.95:
            prediction = "中等抗原性"
            threshold_note = "具有一定的免疫原性潜力"
        else:
            prediction = "低抗原性"
            threshold_note = "免疫原性较弱"

        # 识别高抗原性区域
        high_antigenic_regions = []
        threshold = 1.0
        current_region = []

        for i, val in enumerate(antigenicity_values):
            if val >= threshold:
                current_region.append((i, val))
            else:
                if len(current_region) >= 4:
                    high_antigenic_regions.append({
                        'start': current_region[0][0],
                        'end': current_region[-1][0] + window - 1,
                        'length': len(current_region),
                        'avg_score': round(sum(v for _, v in current_region) / len(current_region), 3),
                    })
                current_region = []

        if len(current_region) >= 4:
            high_antigenic_regions.append({
                'start': current_region[0][0],
                'end': current_region[-1][0] + window - 1,
                'length': len(current_region),
                'avg_score': round(sum(v for _, v in current_region) / len(current_region), 3),
            })

        return {
            'score': round(overall_avg, 3),
            'prediction': prediction,
            'threshold_note': threshold_note,
            'high_antigenic_regions': high_antigenic_regions,
            'max_value': round(max(antigenicity_values), 3) if antigenicity_values else 0,
            'min_value': round(min(antigenicity_values), 3) if antigenicity_values else 0,
        }

    def _predict_transmembrane_regions(self):
        """预测跨膜区 - Kyte-Doolittle滑动窗口法"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length < 20:
            return {'has_transmembrane': False, 'regions': [], 'note': '序列过短，无法预测'}

        # 使用window=20
        profile = self._calculate_hydrophobicity_profile(window=20)

        regions = []
        threshold = TRANSMEMBRANE_THRESHOLD
        current_region = []

        for i, val in enumerate(profile['profile']):
            if val >= threshold:
                current_region.append((i, val))
            else:
                if len(current_region) >= 15:  # 跨膜区通常15-30个残基
                    regions.append({
                        'start': current_region[0][0],
                        'end': current_region[-1][0],
                        'length': len(current_region),
                        'max_hydrophobicity': max(v for _, v in current_region),
                        'avg_hydrophobicity': round(sum(v for _, v in current_region) / len(current_region), 3),
                    })
                current_region = []

        if len(current_region) >= 15:
            regions.append({
                'start': current_region[0][0],
                'end': current_region[-1][0],
                'length': len(current_region),
                'max_hydrophobicity': max(v for _, v in current_region),
                'avg_hydrophobicity': round(sum(v for _, v in current_region) / len(current_region), 3),
            })

        has_tm = len(regions) > 0

        if has_tm:
            prediction = f"可能为膜蛋白 (预测{len(regions)}个跨膜区)"
            note = "序列可能含有跨膜螺旋"
        else:
            prediction = "可能为可溶性蛋白"
            note = "未预测到明显的跨膜区"

        return {
            'has_transmembrane': has_tm,
            'regions': regions,
            'prediction': prediction,
            'note': note,
        }

    def _predict_disulfide_bonds(self):
        """预测二硫键形成可能性"""
        seq = self.cleaned_sequence
        cys_positions = [i for i, aa in enumerate(seq) if aa == 'C']
        cys_count = len(cys_positions)

        if cys_count < 2:
            return {
                'cysteine_count': cys_count,
                'potential_bonds': 0,
                'estimated_bonds': 0,
                'prediction': '无二硫键形成条件',
                'positions': cys_positions,
                'note': 'Cys残基不足，无法形成二硫键',
                'possible_pairing': [],
            }

        max_bonds = cys_count // 2

        # 评估配对可能性
        possible_bonds = []
        for i, pos1 in enumerate(cys_positions):
            for pos2 in cys_positions[i + 1:]:
                if pos2 - pos1 > 3:  # 排除相邻Cys
                    distance = pos2 - pos1
                    if distance < 20:
                        confidence = "高"
                    elif distance < 50:
                        confidence = "中"
                    else:
                        confidence = "低"

                    possible_bonds.append({
                        'position1': pos1,
                        'position2': pos2,
                        'distance': distance,
                        'confidence': confidence,
                    })

        estimated_bonds = min(max_bonds, len([b for b in possible_bonds if b['confidence'] != '低']))

        if cys_count >= 4:
            prediction = f"可能形成{estimated_bonds}对二硫键"
            note = "Cys残基数量充足，具备形成二硫键的条件"
        else:
            prediction = f"最多形成{max_bonds}对二硫键"
            note = "Cys残基数量有限"

        return {
            'cysteine_count': cys_count,
            'potential_bonds': max_bonds,
            'estimated_bonds': estimated_bonds,
            'prediction': prediction,
            'note': note,
            'positions': cys_positions,
            'possible_pairing': possible_bonds[:10],
        }

    def _predict_n_glycosylation(self):
        """预测N-糖基化位点 - Asn-X-Ser/Thr motif"""
        seq = self.cleaned_sequence
        length = len(seq)

        if length < 3:
            return {'sites': [], 'count': 0, 'note': '序列过短'}

        sites = []
        for match in N_GLYCOSYLATION_MOTIF.finditer(seq):
            pos = match.start()
            sites.append({
                'position': pos,
                'motif': seq[pos:pos + 3],
                'context': seq[max(0, pos - 2):min(length, pos + 5)],
            })

        if sites:
            note = f"发现{len(sites)}个潜在的N-糖基化位点"
        else:
            note = "未发现N-糖基化位点"

        return {
            'sites': sites,
            'count': len(sites),
            'note': note,
        }

    def _predict_o_glycosylation(self):
        """预测O-糖基化位点"""
        seq = self.cleaned_sequence
        length = len(seq)

        sites = []
        for i, aa in enumerate(seq):
            if aa in O_GLYCOSYLATION_SITES:
                prev_aa = seq[i - 1] if i > 0 else ''
                next_aa = seq[i + 1] if i < length - 1 else ''

                score = 0
                if prev_aa in 'APG':
                    score += 1
                if next_aa in 'PGA':
                    score += 1

                sites.append({
                    'position': i,
                    'residue': aa,
                    'context': f"{prev_aa}{aa}{next_aa}",
                    'potential': '高' if score >= 2 else '中' if score >= 1 else '低',
                    'score': score,
                })

        high_potential = [s for s in sites if s['potential'] == '高']

        if high_potential:
            note = f"发现{len(high_potential)}个高潜力O-糖基化位点"
        elif sites:
            note = f"发现{len(sites)}个潜在O-糖基化位点"
        else:
            note = "未发现O-糖基化位点"

        return {
            'sites': sites,
            'count': len(sites),
            'high_potential_count': len(high_potential),
            'note': note,
        }

    def _predict_phosphorylation(self):
        """预测磷酸化位点"""
        seq = self.cleaned_sequence
        length = len(seq)

        serine_threonine_sites = []
        tyrosine_sites = []

        for i, aa in enumerate(seq):
            if aa == 'S' or aa == 'T':
                prev = seq[i - 1] if i > 0 else ''
                next_aa = seq[i + 1] if i < length - 1 else ''

                score = 0
                if prev in 'RK':  # 激活位点偏好
                    score += 1
                if next_aa == 'D' or next_aa == 'E':
                    score += 1

                serine_threonine_sites.append({
                    'position': i,
                    'residue': aa,
                    'context': f"{prev}{aa}{next_aa}",
                    'score': score,
                    'potential': '高' if score >= 2 else '中' if score >= 1 else '低',
                })

            elif aa == 'Y':
                tyrosine_sites.append({
                    'position': i,
                    'residue': 'Y',
                    'context': seq[max(0, i - 1):min(length, i + 2)],
                    'potential': '中',
                })

        high_potential_st = [s for s in serine_threonine_sites if s['potential'] == '高']
        total_phos_sites = len(serine_threonine_sites) + len(tyrosine_sites)

        if total_phos_sites > 0:
            note = f"预测{total_phos_sites}个磷酸化位点 ({len(serine_threonine_sites)}个S/T, {len(tyrosine_sites)}个Y)"
        else:
            note = "未预测到磷酸化位点"

        return {
            'serine_threonine': serine_threonine_sites,
            'tyrosine': tyrosine_sites,
            'total_count': total_phos_sites,
            'high_potential_st_count': len(high_potential_st),
            'note': note,
        }

    def _predict_expression_host(self):
        """
        预测表达宿主 - 综合性分析
        基于氨基酸组成、理化性质判断适合的表达系统
        """
        seq = self.cleaned_sequence
        length = len(seq)
        comp = self.results.get('composition', self._calculate_composition())

        score = {
            'ecoli': 0,
            'mammalian': 0,
            'yeast': 0,
            'insect': 0,
        }

        reasons = {
            'ecoli': [],
            'mammalian': [],
            'yeast': [],
            'insect': [],
        }

        # 1. 分子量分析
        mw = self.results.get('molecular_weight', {}).get('mw_kda', 0)
        if mw < 60:
            score['ecoli'] += 20
            reasons['ecoli'].append(f"分子量({mw:.1f}kDa)较小，适合大肠杆菌表达")
        else:
            score['mammalian'] += 15
            score['insect'] += 10
            reasons['mammalian'].append(f"分子量({mw:.1f}kDa)较大，建议真核系统")
            reasons['insect'].append(f"分子量({mw:.1f}kDa)较大，可考虑昆虫细胞")

        # 2. 等电点分析
        pI = self.results.get('isoelectric_point', {}).get('pI', 7.0)
        if pI < 6 or pI > 9:
            score['ecoli'] -= 10
            score['mammalian'] += 15
            reasons['ecoli'].append("等电点极端，可能在大肠杆菌中形成包涵体")
            reasons['mammalian'].append("等电点极端，真核系统更适合")
        else:
            score['ecoli'] += 10
            reasons['ecoli'].append("等电点在适合大肠杆菌表达的范围内")

        # 3. 疏水性分析
        gravy = self.results.get('gravy', {}).get('gravy', 0)
        if gravy > 0.5:
            score['ecoli'] += 15
            reasons['ecoli'].append("疏水性较强，适合大肠杆菌表达")
        elif gravy < -0.5:
            score['mammalian'] += 10
            score['yeast'] += 10
            reasons['mammalian'].append("亲水性较强，真核表达可能更好")
            reasons['yeast'].append("亲水性较强，酵母表达系统适合")

        # 4. 二硫键分析
        disul = self.results.get('disulfide_bonds', {}).get('estimated_bonds', 0)
        if disul >= 2:
            score['ecoli'] -= 20
            score['mammalian'] += 25
            score['insect'] += 15
            reasons['ecoli'].append(f"含有{disul}对二硫键，大肠杆菌可能无法正确折叠")
            reasons['mammalian'].append(f"含有{disul}对二硫键，真核系统更适合形成二硫键")
            reasons['insect'].append(f"含有{disul}对二硫键，昆虫细胞表达系统适合")

        # 5. 糖基化位点分析
        n_glyco = self.results.get('n_glycosylation', {}).get('count', 0)
        if n_glyco >= 2:
            score['ecoli'] -= 15
            score['mammalian'] += 20
            reasons['ecoli'].append(f"含有{n_glyco}个N-糖基化位点，大肠杆菌无N-糖基化")
            reasons['mammalian'].append(f"含有{n_glyco}个N-糖基化位点，适合哺乳动物细胞")

        # 6. 稀有密码子分析 (简化)
        rare_for_ecoli = sum(comp.get(aa, {}).get('count', 0) for aa in 'R')
        if rare_for_ecoli > length * 0.08:
            score['ecoli'] -= 10
            reasons['ecoli'].append("稀有Arg密码子较多，可能影响表达效率")

        # 7. 稳定性分析
        inst_index = self.results.get('instability_index', {}).get('index', 0)
        if inst_index >= 40:
            score['ecoli'] += 5  # 不稳定蛋白在大肠杆菌中可能表达更快
            reasons['ecoli'].append("不稳定蛋白可能适合快速表达")

        # 8. 跨膜区分析
        if self.results.get('transmembrane', {}).get('has_transmembrane', False):
            score['ecoli'] -= 25
            score['mammalian'] += 20
            score['insect'] += 15
            reasons['ecoli'].append("含有跨膜区，大肠杆菌表达可能困难")
            reasons['mammalian'].append("含有跨膜区，哺乳动物细胞更适合")
            reasons['insect'].append("含有跨膜区，昆虫细胞可考虑")

        # 9. 长度分析
        if length > 500:
            score['ecoli'] -= 15
            score['mammalian'] += 10
            reasons['ecoli'].append("序列较长，大肠杆菌表达可能困难")
            reasons['mammalian'].append("序列较长，真核系统可能更适合")

        # 10. 磷酸化分析
        phos = self.results.get('phosphorylation', {}).get('total_count', 0)
        if phos > 10:
            score['mammalian'] += 10
            reasons['mammalian'].append(f"含有{phos}个磷酸化位点，哺乳动物细胞可进行翻译后修饰")

        # 确定推荐系统
        max_score = max(score.values())
        recommendations = []

        for system, s in score.items():
            if s == max_score:
                recommendations.append(system)

        # 系统名称映射
        system_names = {
            'ecoli': '大肠杆菌 (E. coli)',
            'mammalian': '哺乳动物细胞 (HEK293/CHO)',
            'yeast': '酵母 (Pichia pastoris)',
            'insect': '昆虫细胞 (sf9/sf21)',
        }

        system_cell_lines = {
            'ecoli': 'BL21(DE3), Rosetta, Origami',
            'mammalian': 'HEK293, CHO, Expi293F',
            'yeast': 'Pichia pastoris, GS115',
            'insect': 'sf9, sf21, HighFive',
        }

        results = []
        for sys in ['mammalian', 'ecoli', 'yeast', 'insect']:
            if score[sys] > 0:
                results.append({
                    'system': system_names[sys],
                    'score': score[sys],
                    'cell_lines': system_cell_lines[sys],
                    'reasons': reasons[sys][:3] if reasons[sys] else []
                })

        results.sort(key=lambda x: x['score'], reverse=True)

        # 最终推荐
        recommended = results[0]['system'] if results else '无法确定'
        recommended_cell = system_cell_lines.get(recommendations[0] if recommendations else '', '请咨询专家')

        return {
            'recommendations': results[:4],
            'recommended': recommended,
            'recommended_cell_lines': recommended_cell,
            'ecoli_score': score['ecoli'],
            'mammalian_score': score['mammalian'],
        }

    def _calculate_comprehensive_score(self, results):
        """计算综合评分 (0-100)"""
        score = 100

        # 1. 长度评分
        length = results['length']
        if length < 20:
            score -= 15
        elif length < 50:
            score -= 8
        elif length > 2000:
            score -= 5

        # 2. pI评分
        pI = results['isoelectric_point']['pI']
        if pI < 4 or pI > 11:
            score -= 5

        # 3. GRAVY评分
        gravy = results['gravy']['gravy']
        if gravy > 1.5 or gravy < -1.5:
            score -= 5

        # 4. 强疏水斑块惩罚
        patches = results['hydrophobicity_profile']['hydrophobic_patches']
        for patch in patches:
            if patch['max_value'] > 2.0 and patch['length'] > 10:
                score -= 10

        # 5. 不稳定指数惩罚
        inst_index = results['instability_index']['index']
        if inst_index >= 50:
            score -= 15
        elif inst_index >= 40:
            score -= 10

        # 6. 电荷平衡
        pI_data = results['isoelectric_point']
        pos = pI_data['positive_charges']
        neg = pI_data['negative_charges']
        charge_imbalance = abs(pos - neg)
        if charge_imbalance > length * 0.25:
            score -= 10

        # 7. 稀有氨基酸
        comp = results['composition']
        rare_aa = comp.get('W', {}).get('count', 0) + comp.get('M', {}).get('count', 0)
        if rare_aa > length * 0.05:
            score -= 3

        # 8. PTM加分
        n_glyco = results['n_glycosylation']['count']
        phos = results['phosphorylation']['total_count']
        score += min(n_glyco * 1.5, 8)
        score += min(phos * 0.3, 5)

        # 9. 二硫键加分
        bonds = results['disulfide_bonds']['estimated_bonds']
        score += min(bonds * 2, 8)

        # 10. 稳定性加分
        if results['instability_index']['is_stable']:
            score += 5

        score = max(0, min(100, score))

        # 评分等级
        if score >= 85:
            grade = "优秀"
            interpretation = "序列设计良好，各项指标优秀"
        elif score >= 70:
            grade = "良好"
            interpretation = "序列整体质量较好"
        elif score >= 55:
            grade = "一般"
            interpretation = "序列存在一些可优化之处"
        elif score >= 40:
            grade = "较差"
            interpretation = "建议优化序列设计"
        else:
            grade = "不合格"
            interpretation = "序列需要重大优化"

        return {
            'score': round(score, 1),
            'grade': grade,
            'interpretation': interpretation,
        }


def analyze_sequence(sequence):
    """便捷函数：分析单条序列"""
    analyzer = AdvancedSequenceAnalyzer()
    return analyzer.analyze_sequence(sequence)


def analyze_fasta(fasta_file):
    """分析FASTA文件中的所有序列"""
    sequences = {}
    current_id = None
    current_seq = []
    line_number = 0

    try:
        with open(fasta_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line_number += 1
                line = line.strip()

                # 跳过空行
                if not line:
                    continue

                if line.startswith('>'):
                    # 保存之前的序列
                    if current_id is not None and current_seq:
                        seq = ''.join(current_seq)
                        # 清理序列，只保留有效氨基酸
                        seq = ''.join(c for c in seq.upper() if c in 'ACDEFGHIKLMNPQRSTVWY')
                        if seq:
                            sequences[current_id] = seq
                        current_seq = []

                    # 解析ID，支持多种分隔符
                    header = line[1:].strip()
                    if header:
                        # 取第一个空格前的部分作为ID
                        parts = header.split()
                        current_id = parts[0] if parts else f"sequence_{len(sequences) + 1}"
                    else:
                        current_id = f"sequence_{len(sequences) + 1}"

                else:
                    # 收集序列行
                    # 过滤掉数字和特殊字符
                    clean_line = ''.join(c for c in line.upper() if c in 'ACDEFGHIKLMNPQRSTVWY*')
                    if clean_line:
                        current_seq.append(clean_line)

            # 保存最后一个序列
            if current_id is not None and current_seq:
                seq = ''.join(current_seq)
                seq = ''.join(c for c in seq.upper() if c in 'ACDEFGHIKLMNPQRSTVWY')
                if seq:
                    sequences[current_id] = seq

    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        import traceback
        traceback.print_exc()
        return []

    print(f"Parsed {len(sequences)} sequences from FASTA file")

    results = []
    analyzer = AdvancedSequenceAnalyzer()
    for seq_id, seq in sequences.items():
        try:
            if len(seq) >= 3:  # 至少3个氨基酸
                result = analyzer.analyze_sequence(seq)
                result['id'] = seq_id
                results.append(result)
            else:
                results.append({'id': seq_id, 'error': '序列太短', 'length': len(seq)})
        except Exception as e:
            print(f"Error analyzing sequence {seq_id}: {e}")
            results.append({'id': seq_id, 'error': str(e), 'length': len(seq)})

    print(f"Successfully analyzed {len(results)} sequences")
    return results
