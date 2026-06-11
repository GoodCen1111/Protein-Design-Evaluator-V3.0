# -*- coding: utf-8 -*-
"""
PDB Structure Analyzer v3.0
PDB结构分析器

分析蛋白质结构质量，几何参数，以及序列性质
"""

import os
import re
import math
import glob
import numpy as np
from collections import defaultdict
from pathlib import Path


STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

SMALL_MOLECULES = {
    'NA', 'MG', 'K', 'CA', 'FE', 'ZN', 'CU', 'MN', 'CO', 'NI',
    'NAD', 'NADP', 'FAD', 'FMN', 'HEM', 'HEME', 'ATP', 'ADP', 'AMP',
    'HOH', 'WAT', 'H2O', 'SO4', 'PO4', 'GLC', 'EDO', 'PEG',
}


class PDBStructureAnalyzer:
    """PDB结构分析器 v3.0"""

    def __init__(self):
        self.results = []

    def parse_pdb_atoms(self, pdb_path):
        """解析PDB文件中的原子信息"""
        atoms = []
        chains_info = {}

        try:
            with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if not (line.startswith("ATOM") or line.startswith("HETATM")):
                        continue

                    try:
                        atom = {
                            "atom_name": line[12:16].strip(),
                            "res_name": line[17:20].strip(),
                            "chain_id": line[21].strip(),
                            "res_seq": int(line[22:26].strip()),
                            "x": float(line[30:38]),
                            "y": float(line[38:46]),
                            "z": float(line[46:54]),
                            "is_hetatm": line.startswith("HETATM")
                        }

                        try:
                            atom["b_factor"] = float(line[60:66]) if len(line) > 66 else 1.0
                        except:
                            atom["b_factor"] = 1.0

                        atoms.append(atom)

                        if atom["chain_id"] not in chains_info:
                            chains_info[atom["chain_id"]] = {
                                'res_names': set(),
                                'residue_numbers': set()
                            }
                        chains_info[atom["chain_id"]]['res_names'].add(atom["res_name"])
                        chains_info[atom["chain_id"]]['residue_numbers'].add(atom["res_seq"])

                    except:
                        continue
        except:
            pass

        return atoms, chains_info

    def classify_chain(self, chain_id, chains_info, atoms):
        """分类链类型"""
        if chain_id not in chains_info:
            return 'unknown'

        res_names = chains_info[chain_id]['res_names']

        if res_names <= {'HOH', 'WAT', 'H2O', 'DOD'}:
            return 'water'

        protein_res = res_names & STANDARD_AMINO_ACIDS
        if len(protein_res) >= 1:
            return 'protein'

        if res_names & {'DA', 'DT', 'DC', 'DG', 'A', 'T', 'C', 'G', 'U'}:
            return 'nucleic_acid'

        for res in res_names:
            if res in SMALL_MOLECULES:
                return 'small_molecule'

        if len(res_names) <= 3 and not protein_res:
            return 'small_molecule'

        return 'unknown'

    def get_protein_chains(self, atoms, chains_info):
        """获取所有蛋白质链"""
        protein_chains = []

        for chain_id in chains_info:
            chain_type = self.classify_chain(chain_id, chains_info, atoms)
            if chain_type != 'protein':
                continue

            chain_atoms = [a for a in atoms if a['chain_id'] == chain_id]
            res_numbers = sorted(chains_info[chain_id]['residue_numbers'])

            sequence = self._extract_sequence(chain_atoms, chain_id)

            protein_chains.append({
                'chain_id': chain_id,
                'residue_numbers': res_numbers,
                'min_res': min(res_numbers) if res_numbers else 0,
                'max_res': max(res_numbers) if res_numbers else 0,
                'full_length': (max(res_numbers) - min(res_numbers) + 1) if res_numbers else 0,
                'ca_count': len(set(a['res_seq'] for a in chain_atoms if a['atom_name'] == 'CA')),
                'sequence': sequence,
                'atoms': chain_atoms,
            })

        return protein_chains

    def _extract_sequence(self, chain_atoms, chain_id):
        """提取氨基酸序列"""
        aa_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XAA': 'X',
            'UNK': 'X', 'MSE': 'M', 'CSO': 'C', 'PTR': 'Y', 'SEP': 'S',
            'TPO': 'T', 'KCX': 'K', 'CSD': 'C', 'MLY': 'K',
        }

        res_atoms = defaultdict(list)
        for atom in chain_atoms:
            if atom['chain_id'] != chain_id:
                continue
            res_seq = atom['res_seq']
            res_atoms[res_seq].append(atom)

        residues = {}
        for res_seq, atoms in res_atoms.items():
            res_names = set(a['res_name'] for a in atoms)
            if res_names & {'HOH', 'WAT', 'H2O', 'DOD'}:
                continue

            ca_atoms = [a for a in atoms if a['atom_name'] == 'CA']
            if ca_atoms:
                res_name = ca_atoms[0]['res_name']
            else:
                n_atoms = [a for a in atoms if a['atom_name'] == 'N']
                if n_atoms:
                    res_name = n_atoms[0]['res_name']
                else:
                    for name in ['C', 'O', 'CA', 'N']:
                        found = [a for a in atoms if a['atom_name'] == name]
                        if found:
                            res_name = found[0]['res_name']
                            break
                    else:
                        continue

            residues[res_seq] = {'res_name': res_name, 'atom_count': len(atoms)}

        sequence = []
        for res_seq in sorted(residues.keys()):
            res_name = residues[res_seq]['res_name']
            aa = aa_map.get(res_name, 'X')
            sequence.append(aa)

        return ''.join(sequence)

    def detect_linker(self, chain_atoms, chain_id):
        """检测Linker区域"""
        res_by_seq = {}
        for a in chain_atoms:
            if a["atom_name"] == "CA":
                res_seq = a["res_seq"]
                if res_seq not in res_by_seq:
                    res_by_seq[res_seq] = {"res_name": a["res_name"], "b_factor": a["b_factor"]}
                elif a["b_factor"] == 0:
                    res_by_seq[res_seq]["b_factor"] = 0

        residues = sorted(res_by_seq.items(), key=lambda x: x[0])
        best_linker = None

        i = 0
        while i < len(residues):
            if residues[i][1]["res_name"] == "GLY":
                j = i
                while j < len(residues) and residues[j][1]["res_name"] == "GLY":
                    j += 1
                run_len = j - i
                if run_len >= 5:
                    start = residues[i][0]
                    end = residues[j - 1][0]
                    if best_linker is None or run_len > best_linker[2]:
                        best_linker = (start, end, run_len)
                i = j if j > i else i + 1
            else:
                i += 1

        if best_linker is None:
            b_zero = [r for r, info in res_by_seq.items() if info["b_factor"] == 0]
            if len(b_zero) >= 5:
                b_zero.sort()
                best_linker = (min(b_zero), max(b_zero), len(b_zero))

        return best_linker

    def distance(self, a, b):
        return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

    def dihedral(self, p1, p2, p3, p4):
        def vsub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
        def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        def norm(a): return math.sqrt(dot(a, a))

        b1, b2, b3 = vsub(p2, p1), vsub(p3, p2), vsub(p4, p3)
        n1, n2 = cross(b1, b2), cross(b2, b3)
        b2n, n1n, n2n = norm(b2), norm(n1), norm(n2)

        if b2n == 0 or n1n == 0 or n2n == 0:
            return None

        n1u, n2u, b2u = (n1[0]/n1n, n1[1]/n1n, n1[2]/n1n), (n2[0]/n2n, n2[1]/n2n, n2[2]/n2n), (b2[0]/b2n, b2[1]/b2n, b2[2]/b2n)
        m1 = cross(n1u, b2u)
        x, y = dot(n1u, n2u), dot(m1, n2u)
        return math.degrees(math.atan2(y, x))

    def build_atom_map(self, atoms, chain_id=None):
        """构建原子坐标映射"""
        atom_map = defaultdict(dict)
        for a in atoms:
            if chain_id and a["chain_id"] != chain_id:
                continue
            atom_map[(a["chain_id"], a["res_seq"])][a["atom_name"]] = (a["x"], a["y"], a["z"])
        return atom_map

    def compute_geometry(self, atoms, chain_id, res_start, res_end):
        """计算几何参数"""
        atom_map = self.build_atom_map(atoms, chain_id)

        peptide_bonds = []
        for r in range(res_start, res_end + 1):
            c = atom_map.get((chain_id, r), {}).get("C")
            n_next = atom_map.get((chain_id, r + 1), {}).get("N")
            if c and n_next:
                peptide_bonds.append((r, r+1, self.distance(c, n_next)))

        ca_steps = []
        for r in range(res_start, res_end):
            ca_i = atom_map.get((chain_id, r), {}).get("CA")
            ca_next = atom_map.get((chain_id, r + 1), {}).get("CA")
            if ca_i and ca_next:
                ca_steps.append((r, r+1, self.distance(ca_i, ca_next)))

        omegas = []
        for i in range(res_start, res_end + 1):
            ca_im1 = atom_map.get((chain_id, i-1), {}).get("CA")
            c_im1 = atom_map.get((chain_id, i-1), {}).get("C")
            n_i = atom_map.get((chain_id, i), {}).get("N")
            ca_i = atom_map.get((chain_id, i), {}).get("CA")
            if all([ca_im1, c_im1, n_i, ca_i]):
                omega = self.dihedral(ca_im1, c_im1, n_i, ca_i)
                if omega is not None:
                    omegas.append((i-1, i, abs(omega)))

        phi_psi = []
        for r in range(res_start - 1, res_end + 1):
            ca_i = atom_map.get((chain_id, r), {}).get("CA")
            c_i = atom_map.get((chain_id, r), {}).get("C")
            n_i1 = atom_map.get((chain_id, r + 1), {}).get("N")
            ca_i1 = atom_map.get((chain_id, r + 1), {}).get("CA")
            n_i = atom_map.get((chain_id, r), {}).get("N")
            c_im1 = atom_map.get((chain_id, r - 1), {}).get("C")
            if None in [ca_i, c_i, n_i1, ca_i1, n_i, c_im1]:
                phi_psi.append((None, None))
            else:
                phi = self.dihedral(c_im1, n_i, ca_i, c_i)
                psi = self.dihedral(n_i, ca_i, c_i, n_i1)
                phi_psi.append((phi, psi))

        return peptide_bonds, ca_steps, omegas, phi_psi

    def detect_clashes(self, atoms, chain_id, res_start, res_end):
        """检测空间冲突"""
        region_atoms = [a for a in atoms if a["chain_id"] == chain_id and res_start <= a["res_seq"] <= res_end]
        other_atoms = [a for a in atoms if a["chain_id"] == chain_id and not (res_start <= a["res_seq"] <= res_end)]
        other_chains = [a for a in atoms if a["chain_id"] != chain_id]

        severe, mild = 0, 0

        for la in region_atoms:
            la_xyz = (la["x"], la["y"], la["z"])
            la_res = la["res_seq"]

            for oa in other_atoms:
                if abs(la_res - oa["res_seq"]) <= 2:
                    continue
                d = self.distance(la_xyz, (oa["x"], oa["y"], oa["z"]))
                if d < 2.0: severe += 1
                elif d < 2.4: mild += 1

            for oa in other_chains:
                if oa["atom_name"] != "CA":
                    continue
                d = self.distance(la_xyz, (oa["x"], oa["y"], oa["z"]))
                if d < 2.0: severe += 1
                elif d < 2.4: mild += 1

        return severe, mild

    def compute_rg(self, atoms, chain_id, res_start, res_end):
        """计算回转半径"""
        ca_coords = [(a["x"], a["y"], a["z"]) for a in atoms
                     if a["chain_id"] == chain_id and a["atom_name"] == "CA"
                     and res_start <= a["res_seq"] <= res_end]
        if not ca_coords:
            return None
        cx, cy, cz = [sum(c[i] for c in ca_coords) / len(ca_coords) for i in range(3)]
        rg2 = sum((c[0]-cx)**2 + (c[1]-cy)**2 + (c[2]-cz)**2 for c in ca_coords) / len(ca_coords)
        return math.sqrt(rg2)

    def analyze_chain(self, chain_info, atoms):
        """分析单条链"""
        chain_id = chain_info['chain_id']
        chain_atoms = chain_info['atoms']
        sequence = chain_info['sequence']

        res_start = chain_info['min_res']
        res_end = chain_info['max_res']
        full_len = res_end - res_start + 1

        # 获取实际有CA原子的残基数
        ca_atoms = [a for a in chain_atoms if a["atom_name"] == "CA" and res_start <= a["res_seq"] <= res_end]
        actual_residue_count = len(ca_atoms)

        linker = self.detect_linker(chain_atoms, chain_id)
        has_linker = linker is not None

        peptide_bonds, ca_steps, omegas, phi_psi = self.compute_geometry(atoms, chain_id, res_start, res_end)
        severe_clash, mild_clash = self.detect_clashes(atoms, chain_id, res_start, res_end)
        rg = self.compute_rg(atoms, chain_id, res_start, res_end)

        print(f"[Chain {chain_id}] res_start={res_start}, res_end={res_end}, full_len={full_len}")
        print(f"[Chain {chain_id}] actual_residue_count={actual_residue_count}, peptide_bonds={len(peptide_bonds)}, ca_steps={len(ca_steps)}, omegas={len(omegas)}")
        print(f"[Chain {chain_id}] severe_clash={severe_clash}, mild_clash={mild_clash}")

        atom_map = self.build_atom_map(atoms, chain_id)
        n_xyz = atom_map.get((chain_id, res_start), {}).get("N")
        c_before = atom_map.get((chain_id, res_start - 1), {}).get("C")
        c_xyz = atom_map.get((chain_id, res_end), {}).get("C")
        n_after = atom_map.get((chain_id, res_end + 1), {}).get("N")
        left_pep = self.distance(c_before, n_xyz) if (c_before and n_xyz) else None
        right_pep = self.distance(c_xyz, n_after) if (c_xyz and n_after) else None

        ca_coords = [(a["x"], a["y"], a["z"]) for a in chain_atoms
                     if a["atom_name"] == "CA" and res_start <= a["res_seq"] <= res_end]
        if len(ca_coords) >= 2:
            end_to_end = self.distance(ca_coords[0], ca_coords[-1])
            contour = 3.8 * (len(ca_coords) - 1)
            ext_ratio = end_to_end / contour if contour > 0 else None
        else:
            end_to_end = ext_ratio = None

        geo_score = self._score_geometry(peptide_bonds, ca_steps, omegas, severe_clash, mild_clash, left_pep, right_pep)

        print(f"[Chain {chain_id}] geo_score={geo_score}, left_pep={left_pep}, right_pep={right_pep}")

        # 对于膜蛋白或截断结构，使用实际残基数进行评分
        # 避免因残基编号不连续导致的分数偏低问题
        effective_length = actual_residue_count if actual_residue_count > 0 else full_len

        return {
            'chain_id': chain_id,
            'full_length': chain_info['full_length'],
            'actual_residues': actual_residue_count,
            'analyzed_start': res_start,
            'analyzed_end': res_end,
            'analyzed_length': full_len,
            'effective_length': effective_length,
            'sequence': sequence,
            'sequence_length': len(sequence),
            'has_linker': has_linker,
            'linker_region': linker,
            'peptide_bonds': peptide_bonds,
            'ca_steps': ca_steps,
            'omegas': omegas,
            'phi_psi': phi_psi,
            'severe_clash': severe_clash,
            'mild_clash': mild_clash,
            'rg': rg,
            'end_to_end': end_to_end,
            'ext_ratio': ext_ratio,
            'left_peptide': left_pep,
            'right_peptide': right_pep,
            'geo_score': geo_score,
            'is_truncated': actual_residue_count < full_len,  # 是否为截断结构
        }

    def _score_geometry(self, peptide_bonds, ca_steps, omegas, severe_clash, mild_clash, left_pep, right_pep):
        """计算几何评分"""
        score = 100.0

        # 肽键评分
        if peptide_bonds:
            bad = sum(1 for b in peptide_bonds if b[2] < 1.15 or b[2] > 2.2)
            penalty = min(25, bad * 2)  # 降低惩罚力度
            score -= penalty
            if max(b[2] for b in peptide_bonds) > 2.2:
                score -= min(10, (max(b[2] for b in peptide_bonds) - 2.2) * 8)

        # CA步长评分
        if ca_steps:
            bad = sum(1 for s in ca_steps if s[2] < 2.8 or s[2] > 4.5)
            score -= min(15, bad * 2)

        # Omega角评分
        if omegas:
            cis = sum(1 for w in omegas if w[2] < 60)
            bad = sum(1 for w in omegas if abs(w[2] - 180) > 60)
            score -= min(10, cis * 2 + bad * 1.5)

        # 空间冲突评分
        score -= min(20, severe_clash * 1.5 + mild_clash * 0.3)

        # 末端肽键 - 对于截断结构，这不应大幅影响评分
        # 如果是N端或C端残基缺失导致的None，不应惩罚
        if left_pep is not None:
            score -= 3 if left_pep < 1.15 or left_pep > 2.2 else 0
        if right_pep is not None:
            score -= 3 if right_pep < 1.15 or right_pep > 2.2 else 0

        return max(0, min(100, score))

    def analyze_pdb(self, pdb_path):
        """分析PDB文件"""
        atoms, chains_info = self.parse_pdb_atoms(pdb_path)
        if not atoms:
            print(f"[PDB分析] 错误: 无法解析PDB文件 {pdb_path}")
            return None

        protein_chains = self.get_protein_chains(atoms, chains_info)
        if not protein_chains:
            print(f"[PDB分析] 错误: 未找到蛋白质链 chains_info={list(chains_info.keys())}")
            return None

        print(f"[PDB分析] 找到 {len(protein_chains)} 个蛋白质链")

        chain_results = []
        total_score = 0
        total_length = 0
        total_actual = 0

        for chain_info in protein_chains:
            result = self.analyze_chain(chain_info, atoms)
            chain_results.append(result)

            effective_len = result.get('effective_length', result['full_length'])
            geo_score = result.get('geo_score', 0)

            print(f"[PDB分析] Chain {chain_info['chain_id']}: geo_score={geo_score}, effective_len={effective_len}, full_len={result['full_length']}")

            total_score += geo_score * effective_len
            total_length += result['full_length']
            total_actual += effective_len

        print(f"[PDB分析] total_score={total_score}, total_actual={total_actual}")

        # 使用实际残基数计算平均分
        overall_score = total_score / total_actual if total_actual > 0 else 0

        print(f"[PDB分析] overall_score={overall_score}")

        other_chains = {}
        for chain_id in chains_info:
            if chain_id not in [c['chain_id'] for c in protein_chains]:
                chain_type = self.classify_chain(chain_id, chains_info, atoms)
                other_chains[chain_id] = {
                    'type': chain_type,
                    'residue_count': len(chains_info[chain_id]['residue_numbers']),
                    'residues': list(chains_info[chain_id]['res_names'])[:5]
                }

        return {
            'pdb_path': pdb_path,
            'basename': os.path.basename(pdb_path),
            'total_chains': len(chains_info),
            'protein_chains': len(protein_chains),
            'other_chains': other_chains,
            'chain_results': chain_results,
            'overall_score': overall_score,
            'total_length': total_length,
            'total_actual_residues': total_actual,
        }

    def analyze_directory(self, pdb_dir):
        """分析目录中所有PDB文件"""
        results = []
        pdb_files = glob.glob(os.path.join(pdb_dir, "**", "*.pdb"), recursive=True)

        for pdb_file in sorted(pdb_files):
            result = self.analyze_pdb(pdb_file)
            if result:
                results.append(result)

        return sorted(results, key=lambda x: x['overall_score'], reverse=True)

    def compute_all_metrics(self, results):
        """计算所有指标"""
        for r in results:
            for chain in r.get('chain_results', []):
                if chain.get('peptide_bonds'):
                    lengths = [b[2] for b in chain['peptide_bonds']]
                    chain['mean_peptide'] = np.mean(lengths)
                    chain['frac_ok_peptide'] = sum(1 for d in lengths if 1.15 <= d <= 2.2) / len(lengths)
                if chain.get('ca_steps'):
                    steps = [s[2] for s in chain['ca_steps']]
                    chain['mean_ca_step'] = np.mean(steps)
                    chain['ca_step_std'] = np.std(steps)
                    chain['n_bad_ca'] = sum(1 for s in steps if s < 2.8 or s > 4.5)
                if chain.get('omegas'):
                    chain['n_cis_like'] = sum(1 for w in chain['omegas'] if w[2] < 60)
        return results
