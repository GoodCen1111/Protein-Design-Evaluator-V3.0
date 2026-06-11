# -*- coding: utf-8 -*-
"""
ProteinDesignEvaluator v3.0 - 蛋白质设计综合评估工具
主程序入口

功能:
1. 等电点 (pI) 计算 - 影响纯化、电泳、溶解度、表达
2. 亲疏水性分布图 / 平均疏水值 - 检测强疏水斑块
3. 抗原性 / 免疫原性倾向 - 药物蛋白、抗体、疫苗必备
4. 柔性 / 刚性分布 - 判断linker、活性区合理性
5. 不稳定指数 - 预测蛋白在溶液中是否容易降解
6. 脂肪族指数 - 衡量疏水核心、热稳定性
7. 跨膜区预测 - 膜蛋白识别
8. 二硫键可能性预测 - 结构稳定性关键指标
9. N-糖基化位点 - 表达、折叠、稳定性影响
10. O-糖基化位点 - 翻译后修饰
11. 磷酸化位点 - 活性、调控、稳定性
"""

import sys
import os
import traceback
import datetime
import json
import platform
from pathlib import Path

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QPushButton, QLabel, QFileDialog, QMessageBox,
    QTextEdit, QProgressBar, QGroupBox, QCheckBox, QSpinBox,
    QComboBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QScrollArea, QSplitter, QFrame, QStyleFactory, QProgressDialog,
    QDialog, QGridLayout, QTextBrowser, QButtonGroup, QRadioButton,
    QSplashScreen
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QSize, QTimer
from PyQt5.QtGui import QFont, QIcon, QPalette, QColor, QPixmap, QPainter

from src.analyzer import AdvancedSequenceAnalyzer
from src.pdb_analyzer import PDBStructureAnalyzer
from src.visualizer import ProteinVisualizer


class WorkerThread(QThread):
    """后台工作线程"""
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(bool, str)
    result_ready = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)

    def __init__(self, task_type, params):
        super().__init__()
        self.task_type = task_type
        self.params = params
        self._error_sent = False

    def run(self):
        try:
            if self.task_type == "sequence_single":
                self._run_sequence_single()
            elif self.task_type == "sequence_batch":
                self._run_sequence_batch()
            elif self.task_type == "pdb_single":
                self._run_pdb_single()
            elif self.task_type == "pdb_batch":
                self._run_pdb_batch()
            elif self.task_type == "comprehensive":
                self._run_comprehensive()
        except Exception as e:
            if not self._error_sent:
                self._error_sent = True
                traceback.print_exc()
                self.error_occurred.emit(str(e))

    def _run_sequence_single(self):
        analyzer = AdvancedSequenceAnalyzer()
        self.progress.emit(30, "正在分析序列...")
        sequence = self.params.get('sequence', '')
        result = analyzer.analyze_sequence(sequence)
        self.progress.emit(100, "分析完成")
        self.result_ready.emit(result)
        self.finished.emit(True, f"分析完成: {len(sequence)} 残基")

    def _run_sequence_batch(self):
        from src.analyzer import analyze_fasta
        self.progress.emit(30, "正在读取FASTA文件...")
        results = analyze_fasta(self.params['fasta_file'])
        self.progress.emit(100, "分析完成")
        self.result_ready.emit({'results': results, 'type': 'sequence_batch'})
        self.finished.emit(True, f"分析了 {len(results)} 条序列")

    def _run_pdb_single(self):
        analyzer = PDBStructureAnalyzer()
        self.progress.emit(30, "正在解析PDB结构...")
        result = analyzer.analyze_pdb(self.params['pdb_path'])
        if result:
            self.progress.emit(100, "分析完成")
            self.result_ready.emit(result)
            self.finished.emit(True, f"分析完成: {result['basename']}")
        else:
            self.finished.emit(False, "分析失败或无蛋白质链")

    def _run_pdb_batch(self):
        analyzer = PDBStructureAnalyzer()
        self.progress.emit(30, "正在分析目录...")
        results = analyzer.analyze_directory(self.params['pdb_dir'])
        results = analyzer.compute_all_metrics(results)
        self.progress.emit(100, "分析完成")
        self.result_ready.emit({'results': results, 'type': 'pdb_batch'})
        self.finished.emit(True, f"分析了 {len(results)} 个文件")

    def _run_comprehensive(self):
        analyzer_seq = AdvancedSequenceAnalyzer()
        analyzer_pdb = PDBStructureAnalyzer()

        pdb_path = self.params.get('pdb_path', '')
        sequence = self.params.get('sequence', '')

        seq_result = None
        pdb_result = None

        # 如果有PDB文件，先分析PDB获取序列
        if pdb_path:
            self.progress.emit(10, "正在解析PDB结构...")
            pdb_result = analyzer_pdb.analyze_pdb(pdb_path)
            if pdb_result:
                # 从PDB结果中提取第一条链的序列
                chain_results = pdb_result.get('chain_results', [])
                if chain_results:
                    pdb_sequence = chain_results[0].get('sequence', '')
                    if pdb_sequence:
                        sequence = pdb_sequence

        # 如果有序列，进行序列分析
        if sequence:
            self.progress.emit(40, "正在分析序列...")
            seq_result = analyzer_seq.analyze_sequence(sequence)
        elif pdb_result and not sequence:
            # 如果没有手动输入序列，但有PDB，提示用户
            pass

        self.progress.emit(100, "分析完成")

        self.result_ready.emit({
            'sequence_result': seq_result,
            'pdb_result': pdb_result
        })
        self.finished.emit(True, "综合分析完成")


def get_output_dir():
    """获取跨平台的输出目录"""
    if platform.system() == "Windows":
        base = os.path.join(os.environ.get('USERPROFILE', os.path.expanduser("~")), "ProteinDesignEvaluator_Results")
    else:
        base = os.path.join(os.path.expanduser("~"), "ProteinDesignEvaluator_Results")
    os.makedirs(base, exist_ok=True)
    return base


def open_folder(path):
    """跨平台打开文件夹"""
    if platform.system() == "Windows":
        os.startfile(path)
    elif platform.system() == "Darwin":
        import subprocess
        subprocess.run(["open", path])
    else:
        import subprocess
        subprocess.run(["xdg-open", path])


class MainWindow(QMainWindow):
    """主窗口"""

    def __init__(self):
        super().__init__()
        self.output_dir = get_output_dir()
        self.current_output_dir = None
        self.results = None
        self.last_result_type = None
        self.visualizer = ProteinVisualizer(self.output_dir)
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("ProteinDesignEvaluator v3.0 - 蛋白质设计综合评估工具")
        self.setMinimumSize(1100, 800)
        self.setStyleSheet("""
            QMainWindow {
                background-color: #0d1117;
            }
            QLabel {
                color: #c9d1d9;
            }
            QGroupBox {
                border: 1px solid #30363d;
                border-radius: 8px;
                margin-top: 10px;
                font-weight: bold;
                color: #58a6ff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
            QPushButton {
                background-color: #238636;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 8px 16px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #2ea043;
            }
            QPushButton:pressed {
                background-color: #238636;
            }
            QTextBrowser {
                background-color: #161b22;
                color: #c9d1d9;
                border: 1px solid #30363d;
                border-radius: 5px;
            }
            QTabWidget::pane {
                border: 1px solid #30363d;
                border-radius: 5px;
                background-color: #0d1117;
            }
            QTabBar::tab {
                background-color: #161b22;
                color: #8b949e;
                padding: 10px 20px;
                border: 1px solid #30363d;
                border-bottom: none;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
            }
            QTabBar::tab:selected {
                background-color: #0d1117;
                color: #58a6ff;
                border-bottom: 2px solid #58a6ff;
            }
            QProgressBar {
                border: 1px solid #30363d;
                border-radius: 5px;
                background-color: #161b22;
                text-align: center;
                color: #c9d1d9;
            }
            QProgressBar::chunk {
                background-color: #238636;
                border-radius: 4px;
            }
            QLineEdit, QTextEdit {
                background-color: #161b22;
                color: #c9d1d9;
                border: 1px solid #30363d;
                border-radius: 5px;
                padding: 5px;
            }
            QComboBox {
                background-color: #161b22;
                color: #c9d1d9;
                border: 1px solid #30363d;
                border-radius: 5px;
                padding: 5px;
            }
            QComboBox::drop-down {
                border: none;
            }
            QScrollBar:vertical {
                background: #161b22;
                width: 12px;
            }
            QScrollBar::handle:vertical {
                background: #30363d;
                border-radius: 6px;
            }
        """)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        title_layout = QHBoxLayout()
        title = QLabel("ProteinDesignEvaluator v3.0")
        title.setFont(QFont("Segoe UI", 18, QFont.Bold))
        title.setStyleSheet("color: #58a6ff; padding: 10px;")
        subtitle = QLabel("高级蛋白质设计综合评估工具")
        subtitle.setFont(QFont("Segoe UI", 10))
        subtitle.setStyleSheet("color: #8b949e; padding-left: 10px;")
        title_layout.addWidget(title)
        title_layout.addWidget(subtitle)
        title_layout.addStretch()
        layout.addLayout(title_layout)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.tabs.addTab(self._create_sequence_tab(), "序列分析")
        self.tabs.addTab(self._create_pdb_tab(), "PDB结构分析")
        self.tabs.addTab(self._create_comprehensive_tab(), "综合分析")
        self.tabs.addTab(self._create_results_tab(), "结果查看")
        self.tabs.addTab(self._create_help_tab(), "使用说明")

        status_layout = QHBoxLayout()
        self.status_bar = QLabel("就绪 | v3.0")
        self.status_bar.setStyleSheet("background: #161b22; color: #8b949e; padding: 5px; border-radius: 3px;")
        status_layout.addWidget(self.status_bar)
        status_layout.addStretch()

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumWidth(300)
        status_layout.addWidget(self.progress_bar)

        layout.addLayout(status_layout)

    def _create_sequence_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setSpacing(15)

        info = QLabel("""
        <b>高级序列分析功能:</b><br>
        等电点 (pI) | 亲疏水性分布 | 抗原性预测 | 柔性/刚性分布<br>
        不稳定指数 | 脂肪族指数 | 跨膜区预测 | 二硫键预测<br>
        N/O-糖基化位点 | 磷酸化位点
        """)
        info.setStyleSheet("color: #58a6ff; padding: 15px; background: #161b22; border-radius: 8px; border: 1px solid #30363d;")
        layout.addWidget(info)

        input_group = QGroupBox("序列输入")
        input_layout = QVBoxLayout()

        id_layout = QHBoxLayout()
        id_layout.addWidget(QLabel("序列名称:"))
        self.seq_id_input = QTextEdit()
        self.seq_id_input.setPlaceholderText("输入序列名称（可选）")
        self.seq_id_input.setMaximumHeight(35)
        id_layout.addWidget(self.seq_id_input, 1)
        input_layout.addLayout(id_layout)

        seq_layout = QHBoxLayout()
        seq_layout.addWidget(QLabel("氨基酸序列:"))
        self.seq_input = QTextEdit()
        self.seq_input.setPlaceholderText("输入氨基酸序列 (如: MVLSPADKTN...)")
        self.seq_input.setMinimumHeight(100)
        seq_layout.addWidget(self.seq_input, 1)
        input_layout.addLayout(seq_layout)

        btn_layout = QHBoxLayout()
        analyze_btn = QPushButton("分析序列")
        analyze_btn.clicked.connect(self._run_sequence_single)
        analyze_btn.setStyleSheet("background-color: #238636; padding: 12px 24px; font-size: 12pt;")

        clear_btn = QPushButton("清空")
        clear_btn.clicked.connect(self._clear_sequence_input)
        clear_btn.setStyleSheet("background-color: #da3633; padding: 12px 24px; font-size: 12pt;")

        btn_layout.addWidget(analyze_btn)
        btn_layout.addWidget(clear_btn)

        fasta_layout = QHBoxLayout()
        self.fasta_label = QLabel("未选择FASTA文件")
        fasta_btn = QPushButton("选择FASTA")
        fasta_btn.clicked.connect(lambda: self._select_file('fasta', self.fasta_label, "FASTA文件 (*.fasta *.fa *.faa *.fna *.faa *.txt)"))
        batch_btn = QPushButton("批量分析FASTA")
        batch_btn.clicked.connect(self._run_sequence_batch)
        fasta_clear_btn = QPushButton("X")
        fasta_clear_btn.setMaximumWidth(30)
        fasta_clear_btn.setStyleSheet("background-color: #da3633; color: white; border-radius: 3px;")
        fasta_clear_btn.clicked.connect(lambda: self._clear_label(self.fasta_label))
        fasta_layout.addWidget(QLabel("FASTA:"))
        fasta_layout.addWidget(self.fasta_label, 1)
        fasta_layout.addWidget(fasta_btn)
        fasta_layout.addWidget(batch_btn)
        fasta_layout.addWidget(fasta_clear_btn)
        input_layout.addLayout(fasta_layout)

        input_layout.addLayout(btn_layout)
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        result_group = QGroupBox("分析结果")
        result_layout = QVBoxLayout()
        self.sequence_result = QTextBrowser()
        self.sequence_result.setMinimumHeight(400)
        result_layout.addWidget(self.sequence_result)
        result_group.setLayout(result_layout)
        layout.addWidget(result_group)

        return widget

    def _create_pdb_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        info = QLabel("""
        <b>PDB结构分析功能:</b><br>
        多链自动识别 | 蛋白质/小分子/核酸链分类<br>
        几何参数计算 (肽键、CA步长、二面角)<br>
        空间冲突检测 | Ramachandran图
        """)
        info.setStyleSheet("color: #58a6ff; padding: 15px; background: #161b22; border-radius: 8px; border: 1px solid #30363d;")
        layout.addWidget(info)

        file_group = QGroupBox("PDB文件选择")
        file_layout = QVBoxLayout()

        row1 = QHBoxLayout()
        self.pdb_single_label = QLabel("未选择文件")
        single_btn = QPushButton("选择PDB文件")
        single_btn.clicked.connect(lambda: self._select_file('pdb_single', self.pdb_single_label, "PDB文件 (*.pdb *.PDB)"))
        clear_btn = QPushButton("X")
        clear_btn.setMaximumWidth(30)
        clear_btn.clicked.connect(lambda: self._clear_label(self.pdb_single_label))
        row1.addWidget(QLabel("单文件:"))
        row1.addWidget(self.pdb_single_label, 1)
        row1.addWidget(single_btn)
        row1.addWidget(clear_btn)
        file_layout.addLayout(row1)

        row2 = QHBoxLayout()
        self.pdb_batch_label = QLabel("未选择目录")
        batch_btn = QPushButton("选择目录")
        batch_btn.clicked.connect(lambda: self._select_dir('pdb_batch', self.pdb_batch_label))
        clear_btn2 = QPushButton("X")
        clear_btn2.setMaximumWidth(30)
        clear_btn2.clicked.connect(lambda: self._clear_label(self.pdb_batch_label))
        row2.addWidget(QLabel("批量目录:"))
        row2.addWidget(self.pdb_batch_label, 1)
        row2.addWidget(batch_btn)
        row2.addWidget(clear_btn2)
        file_layout.addLayout(row2)

        file_group.setLayout(file_layout)
        layout.addWidget(file_group)

        btn_layout = QHBoxLayout()
        run_btn = QPushButton("开始分析")
        run_btn.clicked.connect(self._run_pdb_analysis)
        run_btn.setStyleSheet("background-color: #238636; padding: 12px 24px; font-size: 12pt;")

        ramachandran_btn = QPushButton("Ramachandran图")
        ramachandran_btn.clicked.connect(self._generate_ramachandran)
        ramachandran_btn.setStyleSheet("background-color: #1f6feb; padding: 12px 24px; font-size: 12pt;")

        btn_layout.addWidget(run_btn)
        btn_layout.addWidget(ramachandran_btn)
        layout.addLayout(btn_layout)

        result_group = QGroupBox("分析结果")
        result_layout = QVBoxLayout()
        self.pdb_result = QTextBrowser()
        self.pdb_result.setMinimumHeight(400)
        result_layout.addWidget(self.pdb_result)
        result_group.setLayout(result_layout)
        layout.addWidget(result_group)

        return widget

    def _create_comprehensive_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        info = QLabel("""
        <b>综合分析功能:</b><br>
        同时分析PDB结构和氨基酸序列，生成完整的评估报告和可视化图表<br>
        包括: 亲疏水性分布图、柔性分布图、抗原性分布图、PTM位点图、综合仪表板
        """)
        info.setStyleSheet("color: #58a6ff; padding: 15px; background: #161b22; border-radius: 8px; border: 1px solid #30363d;")
        layout.addWidget(info)

        input_group = QGroupBox("同时提供PDB和序列")
        input_layout = QVBoxLayout()

        row1 = QHBoxLayout()
        self.combo_pdb_label = QLabel("未选择PDB文件")
        combo_pdb_btn = QPushButton("选择PDB")
        combo_pdb_btn.clicked.connect(lambda: self._select_file('combo_pdb', self.combo_pdb_label, "PDB文件 (*.pdb *.PDB)"))
        combo_pdb_clear_btn = QPushButton("X")
        combo_pdb_clear_btn.setMaximumWidth(30)
        combo_pdb_clear_btn.setStyleSheet("background-color: #da3633; color: white; border-radius: 3px;")
        combo_pdb_clear_btn.clicked.connect(lambda: self._clear_label(self.combo_pdb_label))
        row1.addWidget(QLabel("PDB文件:"))
        row1.addWidget(self.combo_pdb_label, 1)
        row1.addWidget(combo_pdb_btn)
        row1.addWidget(combo_pdb_clear_btn)
        input_layout.addLayout(row1)

        row2 = QHBoxLayout()
        self.combo_seq_label = QLabel("未输入序列")
        combo_seq_btn = QPushButton("输入序列")
        combo_seq_btn.clicked.connect(self._show_sequence_input_dialog)
        combo_seq_clear_btn = QPushButton("X")
        combo_seq_clear_btn.setMaximumWidth(30)
        combo_seq_clear_btn.setStyleSheet("background-color: #da3633; color: white; border-radius: 3px;")
        combo_seq_clear_btn.clicked.connect(lambda: self._clear_label(self.combo_seq_label))
        row2.addWidget(QLabel("序列(可选):"))
        row2.addWidget(self.combo_seq_label, 1)
        row2.addWidget(combo_seq_btn)
        row2.addWidget(combo_seq_clear_btn)
        input_layout.addLayout(row2)

        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        # 分析按钮
        analyze_btn = QPushButton("开始综合分析")
        analyze_btn.clicked.connect(self._run_comprehensive_analysis)
        analyze_btn.setStyleSheet("background-color: #238636; padding: 12px 24px; font-size: 12pt; font-weight: bold;")
        layout.addWidget(analyze_btn)

        btn_layout = QHBoxLayout()
        dashboard_btn = QPushButton("生成综合仪表板")
        dashboard_btn.clicked.connect(self._generate_dashboard)
        dashboard_btn.setStyleSheet("background-color: #238636; padding: 12px 24px; font-size: 12pt;")

        plots_btn = QPushButton("生成所有图表")
        plots_btn.clicked.connect(self._generate_all_plots)
        plots_btn.setStyleSheet("background-color: #a371f7; padding: 12px 24px; font-size: 12pt;")

        btn_layout.addWidget(dashboard_btn)
        btn_layout.addWidget(plots_btn)
        layout.addLayout(btn_layout)

        result_group = QGroupBox("综合分析结果")
        result_layout = QVBoxLayout()
        self.comprehensive_result = QTextBrowser()
        self.comprehensive_result.setMinimumHeight(400)
        result_layout.addWidget(self.comprehensive_result)
        result_group.setLayout(result_layout)
        layout.addWidget(result_group)

        return widget

    def _create_results_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        toolbar = QHBoxLayout()
        export_btn = QPushButton("导出结果")
        export_btn.clicked.connect(self._export_results)
        export_btn.setStyleSheet("background-color: #238636; padding: 10px 20px;")

        open_btn = QPushButton("打开文件夹")
        open_btn.clicked.connect(self._open_output_folder)
        open_btn.setStyleSheet("background-color: #1f6feb; padding: 10px 20px;")

        toolbar.addWidget(export_btn)
        toolbar.addWidget(open_btn)
        toolbar.addStretch()
        layout.addLayout(toolbar)

        self.results_browser = QTextBrowser()
        self.results_browser.setMinimumHeight(500)
        layout.addWidget(self.results_browser)

        return widget

    def _create_help_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        help_text = QTextBrowser()
        help_text.setHtml("""
        <h1 style="color: #58a6ff;">ProteinDesignEvaluator v3.0 使用说明</h1>

        <h2 style="color: #8b949e;">序列分析</h2>
        <p>输入氨基酸序列进行分析，获取以下信息:</p>
        <ul>
            <li><b>等电点 (pI)</b>: 蛋白质净电荷为零时的pH值，影响纯化和溶解度</li>
            <li><b>亲疏水性分布</b>: 沿序列的疏水性变化，检测强疏水斑块</li>
            <li><b>抗原性预测</b>: Kolaskar & Tongaonkar方法，评估免疫原性</li>
            <li><b>柔性/刚性分布</b>: Karplus & Schulz方法，识别结构域边界</li>
            <li><b>不稳定指数</b>: Guruprasad方法，预测蛋白稳定性 (II>40表示不稳定)</li>
            <li><b>脂肪族指数</b>: 衡量热稳定性，>100表示高热稳定性</li>
            <li><b>跨膜区预测</b>: Kyte-Doolittle方法，识别膜蛋白</li>
            <li><b>二硫键预测</b>: Cys残基配对分析</li>
            <li><b>N-糖基化位点</b>: Asn-X-Ser/Thr motif</li>
            <li><b>O-糖基化位点</b>: Ser/Thr残基分析</li>
            <li><b>磷酸化位点</b>: S/T/Y残基磷酸化潜力</li>
        </ul>

        <h2 style="color: #8b949e;">PDB结构分析</h2>
        <p>上传PDB文件进行结构质量评估:</p>
        <ul>
            <li>多链自动识别和分类</li>
            <li>几何参数计算 (肽键长度、CA步长、二面角)</li>
            <li>空间冲突检测</li>
            <li>Ramachandran图生成</li>
        </ul>

        <h2 style="color: #8b949e;">注意事项</h2>
        <ul>
            <li>序列分析仅支持20种标准氨基酸</li>
            <li>PDB文件需要包含完整的原子坐标</li>
            <li>分析结果会自动保存到用户目录</li>
        </ul>
        """)
        help_text.setStyleSheet("""
            QTextBrowser {
                background-color: #0d1117;
                color: #c9d1d9;
                border: none;
                padding: 20px;
            }
            h1 { color: #58a6ff; }
            h2 { color: #8b949e; border-bottom: 1px solid #30363d; padding-bottom: 5px; }
            li { margin: 5px 0; }
        """)
        layout.addWidget(help_text)

        return widget

    def _clear_label(self, label):
        label_text = label.text()
        label.setText("未选择")
        # 清除存储的文件路径
        self._selected_files = getattr(self, '_selected_files', {})
        for mode, path in list(self._selected_files.items()):
            if path == label_text:
                del self._selected_files[mode]
                break

    def _clear_sequence_input(self):
        self.seq_id_input.clear()
        self.seq_input.clear()

    def _select_file(self, mode, label, filter="所有文件 (*.*)"):
        path, _ = QFileDialog.getOpenFileName(self, "选择文件", os.path.expanduser("~"), filter)
        if path:
            label.setText(path)
            self._selected_files = getattr(self, '_selected_files', {})
            self._selected_files[mode] = path

    def _select_dir(self, mode, label):
        path = QFileDialog.getExistingDirectory(self, "选择目录", os.path.expanduser("~"))
        if path:
            label.setText(path)
            self._selected_files = getattr(self, '_selected_files', {})
            self._selected_files[mode] = path

    def _show_sequence_input_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("输入序列")
        dialog.setMinimumSize(500, 400)

        layout = QVBoxLayout(dialog)
        label = QLabel("请输入氨基酸序列:")
        layout.addWidget(label)

        text_edit = QTextEdit()
        text_edit.setPlaceholderText("MVLSPADKTN...")
        layout.addWidget(text_edit, 1)

        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("确定")
        ok_btn.clicked.connect(lambda: self._set_combo_sequence(text_edit.toPlainText(), dialog))
        cancel_btn = QPushButton("取消")
        cancel_btn.clicked.connect(dialog.reject)
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

        dialog.exec_()

    def _set_combo_sequence(self, sequence, dialog):
        self.combo_seq_label.setText(sequence[:50] + "..." if len(sequence) > 50 else sequence)
        self._selected_files = getattr(self, '_selected_files', {})
        self._selected_files['combo_seq'] = sequence
        dialog.accept()

    def _run_comprehensive_analysis(self):
        """运行综合分析"""
        selected_files = getattr(self, '_selected_files', {})

        pdb_path = selected_files.get('combo_pdb', '')
        sequence = selected_files.get('combo_seq', '')

        if not pdb_path and not sequence:
            QMessageBox.warning(self, "警告", "请至少选择PDB文件或输入序列")
            return

        if pdb_path and not os.path.exists(pdb_path):
            QMessageBox.warning(self, "警告", "PDB文件不存在")
            return

        # 根据输入情况显示提示
        if pdb_path and not sequence:
            self.status_bar.setText("正在分析: 从PDB提取序列...")
        elif sequence and not pdb_path:
            self.status_bar.setText("正在分析: 序列分析...")
        else:
            self.status_bar.setText("正在分析: 序列和结构...")

        self._start_worker("comprehensive", {
            'pdb_path': pdb_path if pdb_path else '',
            'sequence': sequence if sequence else ''
        })

    def _run_sequence_single(self):
        seq = self.seq_input.toPlainText().strip()
        if not seq:
            QMessageBox.warning(self, "警告", "请输入氨基酸序列")
            return

        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        seq = ''.join(c for c in seq.upper() if c in valid_aa)

        if len(seq) < 3:
            QMessageBox.warning(self, "警告", "序列太短，至少需要3个氨基酸")
            return

        self._start_worker("sequence_single", {'sequence': seq})

    def _run_sequence_batch(self):
        fasta = getattr(self, 'fasta_label', None)
        if fasta:
            fasta = fasta.text()
            if fasta == "未选择FASTA文件" or not os.path.exists(fasta):
                QMessageBox.warning(self, "警告", "请先选择有效的FASTA文件")
                return
            self._start_worker("sequence_batch", {'fasta_file': fasta})

    def _run_pdb_analysis(self):
        single = self.pdb_single_label.text()
        batch = self.pdb_batch_label.text()

        if single != "未选择文件" and os.path.exists(single):
            self._start_worker("pdb_single", {'pdb_path': single})
        elif batch != "未选择目录" and os.path.exists(batch):
            self._start_worker("pdb_batch", {'pdb_dir': batch})
        else:
            QMessageBox.warning(self, "警告", "请先选择有效的PDB文件或目录")

    def _start_worker(self, task_type, params):
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.status_bar.setText("分析中...")

        self.worker = WorkerThread(task_type, params)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.result_ready.connect(self._on_result_ready)
        self.worker.error_occurred.connect(self._on_error)
        self.worker.start()

    def _on_progress(self, value, status):
        self.progress_bar.setValue(value)
        self.status_bar.setText(status)

    def _on_finished(self, success, message):
        self.progress_bar.setVisible(False)
        self.status_bar.setText(message)
        if success:
            self.status_bar.setStyleSheet("background: #161b22; color: #3fb950; padding: 5px; border-radius: 3px;")
        else:
            self.status_bar.setStyleSheet("background: #161b22; color: #f85149; padding: 5px; border-radius: 3px;")

    def _on_result_ready(self, data):
        try:
            if isinstance(data, dict):
                if 'results' in data and data.get('type') == 'sequence_batch':
                    self._display_batch_sequence_results(data['results'])
                elif 'results' in data and data.get('type') == 'pdb_batch':
                    self._display_batch_pdb_results(data['results'])
                elif 'sequence_result' in data:
                    self._display_comprehensive_results(data)
                elif 'chain_results' in data:
                    self._display_pdb_single_result(data)
                else:
                    self._display_sequence_single_result(data)
            else:
                self._display_sequence_single_result(data)
        except Exception as e:
            traceback.print_exc()
            QMessageBox.warning(self, "警告", f"显示结果时出错:\n{e}")

    def _on_error(self, msg):
        try:
            self.progress_bar.setVisible(False)
            QMessageBox.critical(self, "错误", f"分析失败:\n{msg}")
        except:
            pass

    def _display_sequence_single_result(self, result):
        self.results = result

        pI = result.get('isoelectric_point', {})
        hydro = result.get('hydrophobicity', {})
        inst = result.get('instability_index', {})
        ali = result.get('aliphatic_index', {})
        mw = result.get('molecular_weight', {})
        comp = result.get('composition', {})
        tm = result.get('transmembrane', {})
        disul = result.get('disulfide_bonds', {})
        n_glyco = result.get('n_glycosylation', {})
        o_glyco = result.get('o_glycosylation', {})
        phos = result.get('phosphorylation', {})
        antigen = result.get('antigenicity', {})
        flex = result.get('flexibility', {})
        score = result.get('comprehensive_score', {})
        gravy = result.get('gravy', {})
        expr = result.get('expression_analysis', {})

        text = f"""
        <h2 style='color: #58a6ff;'>序列分析结果</h2>

        <table style='color: #c9d1d9; width: 100%; border-collapse: collapse;'>
        <tr style='background: #161b22;'>
            <td style='padding: 10px; border: 1px solid #30363d;'><b>序列长度</b></td>
            <td style='padding: 10px; border: 1px solid #30363d;'>{result['length']} aa</td>
            <td style='padding: 10px; border: 1px solid #30363d;'><b>分子量</b></td>
            <td style='padding: 10px; border: 1px solid #30363d;'>{mw.get('mw_kda', 0):.2f} kDa</td>
        </tr>
        </table>

        <h3 style='color: #58a6ff;'>综合评分</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <span style='font-size: 36pt; color: {"#3fb950" if score.get("score", 0) >= 75 else "#f0883e" if score.get("score", 0) >= 50 else "#f85149"};'>{score.get('score', 0):.1f}</span>
            <span style='color: #8b949e;'> / 100</span>
            <br><span style='color: {"#3fb950" if score.get("score", 0) >= 75 else "#f0883e" if score.get("score", 0) >= 50 else "#f85149"};'>{score.get('grade', 'N/A')}</span>
            <p style='color: #8b949e; margin-top: 10px;'>{score.get('interpretation', '')}</p>
        </div>

        <h3 style='color: #58a6ff;'>推荐表达系统</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0; border: 1px solid #238636;'>
            <p style='color: #3fb950; font-size: 14pt; margin-top: 0;'><b>推荐: {expr.get('recommended', 'N/A')}</b></p>
            <p style='color: #c9d1d9;'>细胞系: <b>{expr.get('recommended_cell_lines', 'N/A')}</b></p>
        """

        # 添加各系统评分
        if expr.get('recommendations'):
            text += "<p style='color: #8b949e; margin-bottom: 5px;'>各系统评分:</p>"
            for rec in expr.get('recommendations', [])[:4]:
                sys_score = rec.get('score', 0)
                sys_name = rec.get('system', '')
                bar_width = min(max(sys_score, 0), 100)
                color = '#3fb950' if bar_width >= 30 else '#f0883e' if bar_width >= 15 else '#f85149'
                text += f"""
                <div style='margin: 3px 0;'>
                    <span style='color: #c9d1d9; width: 120px; display: inline-block;'>{sys_name.split('(')[0].strip()}:</span>
                    <div style='display: inline-block; width: 150px; background: #30363d; border-radius: 3px;'>
                        <div style='width: {bar_width}%; background: {color}; height: 12px; border-radius: 3px;'></div>
                    </div>
                    <span style='color: {color};'>{sys_score}</span>
                </div>
                """

        text += "</div>"

        text += f"""
        <h3 style='color: #58a6ff;'>等电点分析</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <span style='font-size: 24pt; color: #58a6ff;'>pI = {pI.get('pI', 'N/A')}</span>
            <p style='color: #c9d1d9;'>蛋白质类型: <b>{pI.get('protein_type', 'N/A')}</b></p>
            <p style='color: #8b949e;'>正电荷残基 (K,R,H): {pI.get('positive_charges', 0)}</p>
            <p style='color: #8b949e;'>负电荷残基 (D,E): {pI.get('negative_charges', 0)}</p>
        </div>

        <h3 style='color: #58a6ff;'>亲疏水性 (GRAVY)</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>GRAVY值: <b style='color: #58a6ff;'>{gravy.get('gravy', 0):.3f}</b> - {gravy.get('interpretation', 'N/A')}</p>
            <p>平均疏水指数: <b style='color: #58a6ff;'>{hydro.get('average', 0):.3f}</b> ({hydro.get('status', 'N/A')})</p>
            <p>疏水残基比例: {hydro.get('hydrophobic_percentage', 0):.1f}%</p>
            <p>亲水残基比例: {hydro.get('hydrophilic_percentage', 0):.1f}%</p>
        </div>

        <h3 style='color: #58a6ff;'>稳定性分析</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>不稳定指数: <b style='color: {"#3fb950" if inst.get("is_stable", True) else "#f85149"};'>{inst.get('index', 'N/A')}</b> - {"稳定" if inst.get("is_stable", True) else "不稳定"}</p>
            <p>脂肪族指数: <b style='color: #58a6ff;'>{ali.get('index', 'N/A')}</b> - {ali.get('prediction', 'N/A')}</p>
        </div>

        <h3 style='color: #58a6ff;'>抗原性预测</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>抗原性评分: <b style='color: #a371f7;'>{antigen.get('score', 'N/A')}</b></p>
            <p>预测: <b>{antigen.get('prediction', 'N/A')}</b> - {antigen.get('threshold_note', '')}</p>
        </div>

        <h3 style='color: #58a6ff;'>柔性/刚性分布</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>平均柔性: <b>{flex.get('average', 'N/A')}</b></p>
            <p>柔性区数量: {len(flex.get('flexible_regions', []))}</p>
            <p>刚性区数量: {len(flex.get('rigid_regions', []))}</p>
        </div>

        <h3 style='color: #58a6ff;'>跨膜区预测</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>{tm.get('prediction', 'N/A')}</p>
            <p style='color: #8b949e;'>{tm.get('note', '')}</p>
        </div>

        <h3 style='color: #58a6ff;'>二硫键</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>Cys残基数量: {disul.get('cysteine_count', 0)}</p>
            <p>可能形成的二硫键: <b>{disul.get('estimated_bonds', 0)}</b> 对</p>
            <p style='color: #8b949e;'>{disul.get('note', '')}</p>
        </div>

        <h3 style='color: #58a6ff;'>糖基化位点</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>N-糖基化: <b>{n_glyco.get('count', 0)}</b> 个位点</p>
            <p style='color: #8b949e;'>{n_glyco.get('note', '')}</p>
            <p>O-糖基化: <b>{o_glyco.get('count', 0)}</b> 个位点 ({o_glyco.get('high_potential_count', 0)} 个高潜力)</p>
        </div>

        <h3 style='color: #58a6ff;'>磷酸化位点</h3>
        <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0;'>
            <p>Ser/Thr磷酸化: <b>{len(phos.get('serine_threonine', []))}</b> 个位点</p>
            <p>Tyr磷酸化: <b>{len(phos.get('tyrosine', []))}</b> 个位点</p>
            <p style='color: #8b949e;'>{phos.get('note', '')}</p>
        </div>

        <h3 style='color: #58a6ff;'>氨基酸组成</h3>
        <table style='color: #c9d1d9; width: 100%; border-collapse: collapse;'>
        <tr style='background: #161b22;'>
            <th style='padding: 8px; border: 1px solid #30363d;'>分类</th>
            <th style='padding: 8px; border: 1px solid #30363d;'>数量</th>
            <th style='padding: 8px; border: 1px solid #30363d;'>百分比</th>
        </tr>
        """

        classif = comp.get('classification', {})
        for cat, label in [
            ('hydrophobic', '疏水氨基酸'),
            ('hydrophilic', '亲水氨基酸'),
            ('charged_positive', '正电荷'),
            ('charged_negative', '负电荷'),
            ('polar', '极性氨基酸'),
            ('aromatic', '芳香族'),
            ('aliphatic', '脂肪族'),
        ]:
            count = classif.get(cat, 0)
            pct = count / result['length'] * 100 if result['length'] > 0 else 0
            text += f"<tr><td style='padding: 8px; border: 1px solid #30363d;'>{label}</td><td style='padding: 8px; border: 1px solid #30363d;'>{count}</td><td style='padding: 8px; border: 1px solid #30363d;'>{pct:.1f}%</td></tr>"

        text += "</table>"

        self.sequence_result.setHtml(text)
        self.results_browser.setHtml(text)
        self.last_result_type = 'sequence_single'

    def _display_batch_sequence_results(self, results):
        text = f"<h2 style='color: #58a6ff;'>批量序列分析结果 ({len(results)} 条)</h2>"
        text += """
        <table style='color: #c9d1d9; width: 100%; border-collapse: collapse;'>
        <tr style='background: #161b22;'>
            <th style='padding: 10px; border: 1px solid #30363d;'>排名</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>序列ID</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>长度</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>综合评分</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>等电点</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>稳定性</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>分子量(kDa)</th>
        </tr>
        """

        sorted_results = sorted(results, key=lambda x: x.get('comprehensive_score', {}).get('score', 0), reverse=True)

        for i, r in enumerate(sorted_results[:50]):
            seq_id = r.get('id', 'Unknown')[:30]
            length = r.get('length', 0)
            score = r.get('comprehensive_score', {}).get('score', 0)
            pI = r.get('isoelectric_point', {}).get('pI', 'N/A')
            stable = "Y" if r.get('instability_index', {}).get('is_stable', True) else "N"
            mw = r.get('molecular_weight', {}).get('mw_kda', 0)
            score_color = "#3fb950" if score >= 75 else "#f0883e" if score >= 50 else "#f85149"

            text += f"""
            <tr>
                <td style='padding: 8px; border: 1px solid #30363d;'>{i+1}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{seq_id}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{length}</td>
                <td style='padding: 8px; border: 1px solid #30363d; color: {score_color};'>{score:.1f}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{pI}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{stable}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{mw:.2f}</td>
            </tr>
            """

        text += "</table>"
        self.sequence_result.setHtml(text)
        self.results_browser.setHtml(text)
        self.last_result_type = 'sequence_batch'
        self.results = results

    def _display_pdb_single_result(self, result):
        # 检查是否有截断结构
        has_truncated = any(c.get('is_truncated', False) for c in result.get('chain_results', []))

        text = f"""
        <h2 style='color: #58a6ff;'>PDB分析结果: {result['basename']}</h2>
        """

        # 如果有截断结构，添加提示
        if has_truncated:
            text += """
            <div style='background: #2d333b; padding: 10px; border-radius: 5px; margin: 10px 0; border-left: 4px solid #f0883e;'>
                <p style='color: #f0883e; margin: 0;'><b>注意：</b>检测到可能为截断结构（如膜蛋白截去膜区后的结构），
                评分基于实际观测到的残基计算，N/C端缺失导致的末端不连续不影响评分。</p>
            </div>
            """

        text += f"""
        <table style='color: #c9d1d9; width: 100%;'>
        <tr>
            <td><b>总分:</b></td>
            <td style='color: #3fb950; font-size: 16pt;'>{result['overall_score']:.2f}</td>
            <td><b>蛋白质链数:</b></td>
            <td>{result['protein_chains']}</td>
            <td><b>实际残基数:</b></td>
            <td>{result.get('total_actual_residues', result['total_length'])} aa</td>
        </tr>
        </table>

        <h3 style='color: #58a6ff;'>链详情:</h3>
        <table style='color: #c9d1d9; width: 100%; border-collapse: collapse;'>
        <tr style='background: #161b22;'>
            <th style='padding: 10px; border: 1px solid #30363d;'>链ID</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>实际长度</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>编号范围</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>几何评分</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>严重冲突</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>轻微冲突</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>回转半径</th>
        </tr>
        """

        for chain in result.get('chain_results', []):
            rg = f"{chain.get('rg', 0):.2f}A" if chain.get('rg') else "N/A"
            range_str = f"{chain.get('analyzed_start', '?')}-{chain.get('analyzed_end', '?')}"
            actual = chain.get('actual_residues', chain.get('full_length', 0))
            truncated_note = " (截断)" if chain.get('is_truncated', False) else ""
            text += f"""
            <tr>
                <td style='padding: 8px; border: 1px solid #30363d;'><b>Chain {chain['chain_id']}{truncated_note}</b></td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{actual} aa</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{range_str}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{chain.get('geo_score', 0):.1f}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{chain.get('severe_clash', 0)}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{chain.get('mild_clash', 0)}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{rg}</td>
            </tr>
            """

        text += "</table>"

        self.pdb_result.setHtml(text)
        self.results_browser.setHtml(text)
        self.last_result_type = 'pdb_single'
        self.results = result

    def _display_batch_pdb_results(self, results):
        text = f"<h2 style='color: #58a6ff;'>批量PDB分析结果 ({len(results)} 个文件)</h2>"
        text += """
        <table style='color: #c9d1d9; width: 100%; border-collapse: collapse;'>
        <tr style='background: #161b22;'>
            <th style='padding: 10px; border: 1px solid #30363d;'>排名</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>文件名</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>总分</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>链数</th>
            <th style='padding: 10px; border: 1px solid #30363d;'>总长度</th>
        </tr>
        """

        for i, r in enumerate(results[:30]):
            score_color = "#3fb950" if r['overall_score'] >= 80 else "#f0883e" if r['overall_score'] >= 60 else "#f85149"
            text += f"""
            <tr>
                <td style='padding: 8px; border: 1px solid #30363d;'>{i+1}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{r['basename'][:40]}</td>
                <td style='padding: 8px; border: 1px solid #30363d; color: {score_color};'>{r['overall_score']:.2f}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{r['protein_chains']}</td>
                <td style='padding: 8px; border: 1px solid #30363d;'>{r['total_length']}</td>
            </tr>
            """

        text += "</table>"
        self.pdb_result.setHtml(text)
        self.results_browser.setHtml(text)
        self.last_result_type = 'pdb_batch'
        self.results = results

    def _display_comprehensive_results(self, data):
        seq_result = data.get('sequence_result', {})
        pdb_result = data.get('pdb_result')

        text = "<h2 style='color: #58a6ff;'>综合分析报告</h2>"

        has_seq_result = seq_result and isinstance(seq_result, dict) and 'length' in seq_result
        has_pdb_result = pdb_result and isinstance(pdb_result, dict) and 'chain_results' in pdb_result

        if has_seq_result:
            score = seq_result.get('comprehensive_score', {})
            text += f"""
            <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0; border: 1px solid #30363d;'>
                <h3 style='color: #58a6ff; margin-top: 0;'>序列分析结果</h3>
                <p>序列长度: <b>{seq_result.get('length', 0)}</b> aa</p>
                <p>综合评分: <b style='color: #3fb950; font-size: 14pt;'>{score.get('score', 0):.1f}</b> / 100 ({score.get('grade', 'N/A')})</p>
                <p>等电点 (pI): <b>{seq_result.get('isoelectric_point', {}).get('pI', 'N/A')}</b></p>
                <p>分子量: <b>{seq_result.get('molecular_weight', {}).get('mw_kda', 0):.2f} kDa</b></p>
                <p>稳定性: <b style='color: {"#3fb950" if seq_result.get("instability_index", {}).get("is_stable", True) else "#f85149"};'>{"稳定" if seq_result.get("instability_index", {}).get("is_stable", True) else "不稳定"}</b>
                   (不稳定指数: {seq_result.get('instability_index', {}).get('index', 'N/A')})</p>
                <p>脂肪族指数: <b>{seq_result.get('aliphatic_index', {}).get('index', 'N/A')}</b>
                   ({seq_result.get('aliphatic_index', {}).get('prediction', '')})</p>
            </div>
            """

        if has_pdb_result:
            chain_info = pdb_result.get('chain_results', [])[0] if pdb_result.get('chain_results') else {}
            text += f"""
            <div style='background: #161b22; padding: 15px; border-radius: 8px; margin: 10px 0; border: 1px solid #30363d;'>
                <h3 style='color: #58a6ff; margin-top: 0;'>结构分析结果</h3>
                <p>文件名: <b>{pdb_result.get('basename', 'N/A')}</b></p>
                <p>蛋白质链数: <b>{pdb_result.get('protein_chains', 0)}</b></p>
                <p>结构评分: <b style='color: #58a6ff; font-size: 14pt;'>{pdb_result.get('overall_score', 0):.1f}</b> / 100</p>
                <p>实际残基数: <b>{pdb_result.get('total_actual_residues', pdb_result.get('total_length', 0))}</b> aa</p>
            """
            if chain_info:
                rg_val = chain_info.get('rg')
                rg_str = f"{rg_val:.2f} A" if rg_val else "N/A"
                text += f"""
                <p>回转半径: <b>{rg_str}</b></p>
                <p>空间冲突: 严重 <b style='color: #f85149;'>{chain_info.get('severe_clash', 0)}</b>,
                            轻微 <b style='color: #f0883e;'>{chain_info.get('mild_clash', 0)}</b></p>
                """
                # 如果是截断结构，添加提示
                if chain_info.get('is_truncated', False):
                    text += """<p style='color: #f0883e;'><b>注意：</b>可能为截断结构（如膜蛋白），N/C端缺失不影响评分</p>"""
            text += "</div>"

        if not has_seq_result and not has_pdb_result:
            text += "<p style='color: #f85149;'>未找到有效的分析结果</p>"

        self.comprehensive_result.setHtml(text)
        self.results = data

    def _generate_ramachandran(self):
        if not self.results or self.last_result_type != 'pdb_single':
            QMessageBox.warning(self, "警告", "请先进行PDB单文件分析")
            return

        try:
            phi_psi_list = []
            for chain in self.results.get('chain_results', []):
                if chain.get('phi_psi'):
                    phi_psi_list = chain['phi_psi']
                    break

            if phi_psi_list:
                path = self.visualizer.plot_ramachandran(
                    phi_psi_list,
                    title=f"Ramachandran - {self.results.get('basename', 'PDB')}"
                )
                if path:
                    self.current_output_dir = self.output_dir
                    reply = QMessageBox.question(self, "成功", f"Ramachandran图已保存\n{path}\n\n是否打开查看?",
                                              QMessageBox.Yes | QMessageBox.No)
                    if reply == QMessageBox.Yes:
                        open_folder(path)
            else:
                QMessageBox.warning(self, "警告", "未找到phi/psi数据")

        except Exception as e:
            traceback.print_exc()
            QMessageBox.critical(self, "错误", f"生成失败:\n{e}")

    def _generate_dashboard(self):
        if not self.results:
            QMessageBox.warning(self, "警告", "请先进行序列分析")
            return

        try:
            if isinstance(self.results, dict) and 'sequence_result' in self.results:
                seq_result = self.results['sequence_result']
                pdb_result = self.results.get('pdb_result')
            elif isinstance(self.results, dict) and 'length' in self.results:
                seq_result = self.results
                pdb_result = None
            else:
                QMessageBox.warning(self, "警告", "请先进行有效的分析")
                return

            path = self.visualizer.create_analysis_dashboard(seq_result, pdb_result)

            if path:
                self.current_output_dir = self.output_dir
                reply = QMessageBox.question(self, "成功", f"综合仪表板已保存\n{path}\n\n是否打开查看?",
                                          QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.Yes:
                    open_folder(path)

        except Exception as e:
            traceback.print_exc()
            QMessageBox.critical(self, "错误", f"生成失败:\n{e}")

    def _generate_all_plots(self):
        if not self.results:
            QMessageBox.warning(self, "警告", "请先进行序列分析")
            return

        try:
            if isinstance(self.results, dict) and 'sequence_result' in self.results:
                seq_result = self.results['sequence_result']
            elif isinstance(self.results, dict) and 'length' in self.results:
                seq_result = self.results
            else:
                QMessageBox.warning(self, "警告", "请先进行有效的分析")
                return

            saved_files = []

            profile = seq_result.get('hydrophobicity_profile', {})
            if profile.get('profile'):
                path = self.visualizer.plot_hydrophobicity_profile(
                    seq_result['sequence'],
                    profile['profile'],
                    profile.get('hydrophobic_patches', []),
                    title=f"亲疏水性分布 - {seq_result.get('id', 'Sequence')[:20]}"
                )
                if path:
                    saved_files.append(path)

            flex = seq_result.get('flexibility', {})
            if flex.get('profile'):
                path = self.visualizer.plot_flexibility_profile(
                    seq_result['sequence'],
                    flex['profile'],
                    flex.get('flexible_regions', []),
                    flex.get('rigid_regions', []),
                    title="柔性/刚性分布"
                )
                if path:
                    saved_files.append(path)

            sequence = seq_result.get('sequence', '')
            if len(sequence) > 6:
                window = 6
                antigenicity_values = []
                ANTIGENICITY_SCALE = {
                    'A': 1.064, 'R': 0.644, 'N': 0.803, 'D': 0.858, 'C': 1.046,
                    'E': 0.858, 'Q': 0.809, 'G': 0.874, 'H': 0.864, 'I': 1.152,
                    'L': 1.236, 'K': 0.644, 'M': 1.303, 'F': 1.268, 'P': 0.858,
                    'S': 0.941, 'T': 0.811, 'W': 1.089, 'Y': 1.064, 'V': 1.268,
                }
                for i in range(len(sequence) - window + 1):
                    window_seq = sequence[i:i + window]
                    avg = sum(ANTIGENICITY_SCALE.get(aa, 1.0) for aa in window_seq) / window
                    antigenicity_values.append(avg)

                if antigenicity_values:
                    path = self.visualizer.plot_antigenicity_profile(
                        sequence, antigenicity_values,
                        seq_result.get('antigenicity', {}).get('high_antigenic_regions', [])
                    )
                    if path:
                        saved_files.append(path)

            path = self.visualizer.plot_ptm_sites(
                seq_result.get('sequence', ''),
                seq_result.get('phosphorylation'),
                seq_result.get('n_glycosylation'),
                seq_result.get('o_glycosylation')
            )
            if path:
                saved_files.append(path)

            path = self.visualizer.plot_amino_acid_composition(seq_result.get('composition', {}))
            if path:
                saved_files.append(path)

            if saved_files:
                files_list = '\n'.join([os.path.basename(f) for f in saved_files])
                reply = QMessageBox.question(
                    self, "成功",
                    f"已生成 {len(saved_files)} 个图表:\n{files_list}\n\n是否打开图表文件夹?",
                    QMessageBox.Yes | QMessageBox.No
                )
                if reply == QMessageBox.Yes:
                    self._open_output_folder()
            else:
                QMessageBox.warning(self, "警告", "未能生成图表")

        except Exception as e:
            traceback.print_exc()
            QMessageBox.critical(self, "错误", f"生成失败:\n{e}")

    def _export_results(self):
        if self.results is None:
            QMessageBox.warning(self, "警告", "没有可导出的结果")
            return

        # 创建导出对话框
        dialog = QDialog(self)
        dialog.setWindowTitle("导出格式选择")
        dialog.setMinimumSize(400, 150)

        # 设置对话框背景
        dialog.setStyleSheet("""
            QDialog {
                background-color: #0d1117;
            }
            QLabel {
                color: #c9d1d9;
                font-size: 12pt;
            }
        """)

        layout = QVBoxLayout(dialog)

        label = QLabel("请选择导出格式:")
        label.setStyleSheet("color: #c9d1d9; font-size: 12pt; padding: 10px;")
        layout.addWidget(label)

        btn_layout = QHBoxLayout()

        btn_json = QPushButton("JSON")
        btn_json.setStyleSheet("background-color: #238636; color: white; padding: 10px 20px; border-radius: 5px; font-size: 11pt;")
        btn_json.setMinimumHeight(40)

        btn_csv = QPushButton("CSV")
        btn_csv.setStyleSheet("background-color: #1f6feb; color: white; padding: 10px 20px; border-radius: 5px; font-size: 11pt;")
        btn_csv.setMinimumHeight(40)

        btn_txt = QPushButton("TXT")
        btn_txt.setStyleSheet("background-color: #a371f7; color: white; padding: 10px 20px; border-radius: 5px; font-size: 11pt;")
        btn_txt.setMinimumHeight(40)

        btn_cancel = QPushButton("取消")
        btn_cancel.setStyleSheet("background-color: #6e7681; color: white; padding: 10px 20px; border-radius: 5px; font-size: 11pt;")
        btn_cancel.setMinimumHeight(40)

        btn_layout.addWidget(btn_json)
        btn_layout.addWidget(btn_csv)
        btn_layout.addWidget(btn_txt)
        btn_layout.addWidget(btn_cancel)

        layout.addLayout(btn_layout)

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        def do_export(ext):
            filepath = os.path.join(self.output_dir, f"analysis_{timestamp}.{ext}")
            try:
                if ext == 'json':
                    content = json.dumps(self.results, indent=2, ensure_ascii=False)
                elif ext == 'csv':
                    content = self._generate_csv()
                    content = '\ufeff' + content
                else:
                    content = self._generate_txt()

                with open(filepath, 'w', encoding='utf-8-sig') as f:
                    f.write(content)

                dialog.accept()
                QMessageBox.information(self, "成功", f"结果已保存到:\n{filepath}")
                self.current_output_dir = self.output_dir

            except Exception as e:
                dialog.accept()
                QMessageBox.critical(self, "错误", f"保存失败:\n{e}")

        btn_json.clicked.connect(lambda: do_export('json'))
        btn_csv.clicked.connect(lambda: do_export('csv'))
        btn_txt.clicked.connect(lambda: do_export('txt'))
        btn_cancel.clicked.connect(dialog.reject)

        dialog.exec_()

    def _generate_csv(self):
        lines = ["# ProteinDesignEvaluator v3.0 分析结果"]

        if isinstance(self.results, dict) and 'length' in self.results:
            r = self.results
            lines.append(f"# 序列ID,{r.get('id', 'Unknown')}")
            lines.append(f"# 序列长度,{r.get('length', 0)}")
            lines.append(f"# 综合评分,{r.get('comprehensive_score', {}).get('score', 0):.2f}")
            lines.append(f"# 等电点,{r.get('isoelectric_point', {}).get('pI', 'N/A')}")
            lines.append(f"# 分子量(kDa),{r.get('molecular_weight', {}).get('mw_kda', 0):.2f}")
            lines.append(f"# 稳定性,{r.get('instability_index', {}).get('prediction', 'N/A')}")
            lines.append(f"# 脂肪族指数,{r.get('aliphatic_index', {}).get('index', 0):.2f}")
            lines.append(f"# N-糖基化位点,{r.get('n_glycosylation', {}).get('count', 0)}")
            lines.append(f"# 磷酸化位点,{r.get('phosphorylation', {}).get('total_count', 0)}")
            lines.append(f"# 二硫键,{r.get('disulfide_bonds', {}).get('estimated_bonds', 0)}")
            lines.append(f"# 序列,{r.get('sequence', '')}")

        return '\n'.join(lines)

    def _generate_txt(self):
        lines = []
        lines.append("=" * 60)
        lines.append("  ProteinDesignEvaluator v3.0 分析报告")
        lines.append("=" * 60)
        lines.append(f"生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")

        if isinstance(self.results, dict) and 'length' in self.results:
            r = self.results
            lines.append("【序列分析结果】")
            lines.append("-" * 40)
            lines.append(f"序列ID: {r.get('id', 'Unknown')}")
            lines.append(f"序列长度: {r.get('length', 0)} aa")
            lines.append(f"综合评分: {r.get('comprehensive_score', {}).get('score', 0):.1f}")
            lines.append(f"评分等级: {r.get('comprehensive_score', {}).get('grade', 'N/A')}")
            lines.append("")
            lines.append("【关键指标】")
            lines.append(f"等电点: {r.get('isoelectric_point', {}).get('pI', 'N/A')}")
            lines.append(f"分子量: {r.get('molecular_weight', {}).get('mw_kda', 0):.2f} kDa")
            lines.append(f"不稳定指数: {r.get('instability_index', {}).get('index', 'N/A')} - {r.get('instability_index', {}).get('prediction', '')}")
            lines.append(f"脂肪族指数: {r.get('aliphatic_index', {}).get('index', 'N/A')} - {r.get('aliphatic_index', {}).get('prediction', '')}")
            lines.append("")
            lines.append("【翻译后修饰】")
            lines.append(f"N-糖基化: {r.get('n_glycosylation', {}).get('count', 0)} 个位点")
            lines.append(f"O-糖基化: {r.get('o_glycosylation', {}).get('count', 0)} 个位点")
            lines.append(f"磷酸化: {r.get('phosphorylation', {}).get('total_count', 0)} 个位点")
            lines.append(f"二硫键: {r.get('disulfide_bonds', {}).get('estimated_bonds', 0)} 对")

        lines.append("")
        lines.append("=" * 60)
        return '\n'.join(lines)

    def _open_output_folder(self):
        try:
            target = self.current_output_dir if self.current_output_dir else self.output_dir
            if not os.path.exists(target):
                os.makedirs(target, exist_ok=True)
            open_folder(target)
        except Exception as e:
            QMessageBox.warning(self, "错误", f"无法打开: {e}")


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("ProteinDesignEvaluator")
    app.setOrganizationName("ProteinDesignEvaluator")

    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
