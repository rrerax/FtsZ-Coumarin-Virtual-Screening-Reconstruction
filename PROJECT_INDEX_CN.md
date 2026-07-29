# 项目索引

本索引用于快速定位 FtsZ 香豆素类似物虚拟筛选复现项目中的主要文件、结果图和结果表。

## 推荐阅读顺序

| 顺序 | 文件 | 用途 |
|---:|---|---|
| 1 | `README.md` | 项目总览、流程、关键图和复现命令 |
| 2 | `PROJECT_INDEX_CN.md` | 中文路径索引 |
| 3 | `docs/technical_report_cn.md` | 中文技术报告，包含方法、结果和局限性 |
| 4 | `docs/technical_report.md` | 英文技术报告，适合 GitHub 展示和技术审阅 |
| 5 | `notebooks/ftsZ_reconstruction.ipynb` | Notebook 形式的表格和图形总览 |
| 6 | `docs/enopt_extension.md` | EnOpt-style receptor ensemble reranking 说明 |
| 7 | `docs/data_availability.md` | 完整本地复跑所需的原始输入文件结构 |

## 主要结果文件

| 结果 | 路径 |
|---|---|
| 项目汇总图 | `results/figures/workflow_summary.png` |
| BP1/BP2 score 分布图 | `results/figures/score_distribution_by_pocket.png` |
| Top 候选分子结构图 | `results/figures/top_ligand_structures.png` |
| 旧 score 与 Vina 复跑 score 对比图 | `results/figures/vina_rerun_vs_original_scores.png` |
| BP1/BP2 重叠候选分子图 | `results/figures/bp_score_scatter_overlap.png` |
| BP1/BP2 结合模式接触图 | `results/figures/binding_mode_top_hits.png` |
| EnOpt-style reranking 图 | `results/figures/enopt_style_reranking_vs_legacy.png` |
| 构象 docking score heatmap | `results/figures/ensemble_score_matrix_heatmap.png` |
| 构象权重图 | `results/figures/enopt_style_conformation_weights.png` |
| 全量 Vina 复跑结果表 | `results/tables/vina_rerun_unique_structures.csv` |
| 多构象 docking score 表 | `results/tables/ensemble_vina_scores.csv` |
| EnOpt-style Top hit 表 | `results/tables/enopt_style_top_hits.csv` |
| EnOpt-style score matrix | `results/tables/enopt_style_score_matrix.csv` |
| 最优复跑分子的近邻残基表 | `results/tables/binding_contacts_top_hits.csv` |
| EnOpt-ready score matrix | `data/processed/enopt_score_matrix.csv` |

## 文件夹说明

| 文件夹 | 内容 |
|---|---|
| `data/processed/` | 清洗后的候选分子表、pocket box 和 score matrix |
| `docs/` | 中英文技术报告、数据审计和数据可用性说明 |
| `notebooks/` | Notebook 总览 |
| `results/figures/` | README 和报告中使用的结果图 |
| `results/tables/` | 清洗、复跑、相关性和近邻残基分析得到的 CSV 表格 |
| `results/pymol/` | PyMOL-ready 脚本，用于在本地恢复 receptor/pose 后进一步渲染 |
| `scripts/` | 数据清洗、口袋重建、对接准备、Vina 复跑、接触分析和绘图脚本 |

## 方法边界

- 本项目是 FtsZ/coumarin HTVS 流程复现和复跑整理；
- docking score 用于计算筛选和候选优先级判断，不等于实验活性；
- BP1/BP2 是两个不同 binding pockets，不是同一 pocket 的多个构象；
- 当前 EnOpt-style 扩展已经加入 5 个 FtsZ receptor conformations 和 450 个 ensemble docking jobs；
- 当前 reranking 是 conformation-weighted consensus，不是严格 supervised EnOpt 模型；
- 若要完整复跑，需要按 `docs/data_availability.md` 恢复原始 Excel、receptor PDB 和 pocket detection 输出。

## 快速复现命令

只重生公开图表：

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/draw_top_ligands.py
python scripts/analyze_enopt_style.py
python scripts/make_summary_figure.py
```

恢复原始输入后的完整本地复跑：

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
micromamba create -y -p .mamba_vina -f environment-vina.yml

python scripts/clean_candidates.py
python scripts/derive_pocket_boxes.py
python scripts/prepare_docking_inputs.py --top-n 0 --max-heavy-atoms 100
python scripts/run_vina_batch.py --exhaustiveness 4 --cpu 4
python scripts/analyze_binding_contacts.py
python scripts/build_receptor_ensemble.py --target-ca-rmsd 0.65 --modes 2
python scripts/run_ensemble_vina.py --exhaustiveness 3 --cpu 4
python scripts/analyze_enopt_style.py
python scripts/draw_top_ligands.py
python scripts/draw_binding_mode_views.py
python scripts/make_summary_figure.py
```
