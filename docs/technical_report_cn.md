# 技术报告：FtsZ 香豆素类似物 HTVS 流程复现与复跑

## 摘要

本项目将既有 FtsZ 高通量虚拟筛选资料整理为一套可复现的 CADD 分析流程。输入资料包括 BP1、BP2 两个结合口袋的候选分子表、口袋检测输出和部分对接结果；本项目将其转化为结构化数据、可复跑脚本和可审阅图表。

当前工作包括：候选表解析、SMILES 标准化、去盐与重复结构合并、BP1/BP2 口袋盒子重建、90 个 cleaned unique structures AutoDock Vina 全量复跑、旧分数与复跑分数对比、最优复跑分子近邻残基分析，以及基于 5 个 FtsZ receptor conformations 的 EnOpt-style 多构象重打分。

> 方法边界：本项目应表述为 FtsZ HTVS 流程复现、多构象 docking pilot 与 EnOpt-style consensus reranking。由于当前资料没有实验 active/decoy 标签，因此不应解读为严格的 supervised EnOpt 模型或实验活性验证。

## 项目背景

原始研究关注 FtsZ 蛋白的潜在抑制剂筛选，候选小分子主要来自多个公开/商业化分子库，并以 coumarin analogue 相关骨架为中心进行虚拟筛选。本项目目标是将原本分散在论文说明、Excel 结果表、pocket detection 输出和 docking 输出中的信息重建为一个结构化计算流程，便于方法复核、结果追踪和后续扩展。

## 数据整理结果

| 指标 | BP1 | BP2 | 合计 / 说明 |
|---|---:|---:|---|
| 原始方法层面记录 | 101 | 103 | 204 条 |
| 有效 SMILES 记录 | 101 | 103 | 204 / 204 |
| 去重后 compound ID | 55 | 72 | 按 pocket 统计 |
| 原始 docking score 最小值 | -7.396 | -9.116 | kcal/mol，越低越好 |
| 原始 docking score 中位数 | -6.980 | -7.429 | kcal/mol |
| 原始 docking score 最大值 | -6.318 | -5.969 | kcal/mol |
| 去盐 + 结构去重后结构数 |  |  | 90 个 unique structures |
| EnOpt-ready matrix |  |  | 88 行 |
| receptor ensemble |  |  | 5 个 FtsZ 构象 |
| multi-conformation docking |  |  | 450 / 450 完成 |

主要清洗逻辑包括：

- 解析旧 Excel 表格中实际对应的 compound ID、原始 docking score、Tanimoto、SMILES 与 fingerprint method 字段；
- 使用 RDKit 校验并 canonicalize SMILES；
- 去除盐和 counter-ion，按 canonical desalted SMILES 合并重复结构；
- 保留不同数据库来源信息，避免同一结构在 ChEMBL、DrugBank、HMDB、MolPort、ZINC 等来源中重复计数；
- 生成方法层面表、compound 层面表、结构层面表和 compound × pocket score matrix。

## 口袋重建

由于输入资料中未找到原始 Vina config 或 ezSMDock 批量配置，本项目根据 pocket detection 输出和文献/原始说明中提到的关键残基重建 BP1/BP2 docking box。该 box 用于本次复现分析，并标注为 estimated docking box。

| Pocket | 口袋来源 | Center X | Center Y | Center Z | Box size | 残基锚点 |
|---|---|---:|---:|---:|---:|---|
| BP1 | `ezPocket_fpocket3.tsv / Pocket1` | -50.6447 | 35.0742 | 15.2025 | 22.5 Å × 22.5 Å × 22.5 Å | Ala39B, Asn41B, Leu47B, Met49B, Ser50B, Lys55B |
| BP2 | `ezPocket_fconv.tsv / Pocket0` | -63.4691 | 20.8178 | 1.7349 | 22.5 Å × 22.5 Å × 22.5 Å | Met169B, Glu185B, Asn189B, Ile225B, Gly226B, Ser227B, Arg304B |

## Vina 复跑结果

对清洗、去盐和结构去重后的 90 个 unique structures 进行 AutoDock Vina 全量复跑。受体使用 FtsZ receptor PDB；配体由 SMILES 重新生成 3D 构象，并转换为 PDBQT。Vina / OpenBabel 环境通过 `environment-vina.yml` 独立记录。其中一个 BP2 大分子使用 loose RDKit 3D-coordinate fallback 完成构象生成，全量复跑最终完成。

| Pocket   |   复跑数量 |   成功数量 | 最优复跑 compound   | 原始 score   | 复跑 score   | 分数变化   |
|:---------|-----------:|-----------:|:--------------------|:-------------|:-------------|:-----------|
| BP1      |         37 |         37 | DB00776             | -6.883       | -8.02        | -1.137     |
| BP2      |         53 |         53 | 686393;CHEMBL392451 | -8.211       | -10.76       | -2.549     |
| 合计     |         90 |         90 |                     |              |              |            |

分数相关性结果：

| 范围    |   n |   Pearson r |   Pearson p |   Spearman r |   Spearman p |
|:--------|----:|------------:|------------:|-------------:|-------------:|
| BP1     |  37 |       0.528 |      0.0008 |        0.543 |       0.0005 |
| BP2     |  53 |       0.743 |      0      |        0.771 |       0      |
| overall |  90 |       0.737 |      0      |        0.737 |       0      |

解释建议：

- Overall 相关性较高，说明重建流程与旧表格分数整体方向一致；
- BP2 的 Pearson 相关性更稳定，Top hit 在复跑中仍保持较强 docking score；
- BP1 的复跑排序变化较明显，可能与 docking box 重建、配体构象、质子化/电荷处理、随机种子或原始参数缺失有关；
- 该结果适合作为复现一致性与结构解释分析，不应作为实验活性证据。

## 结合模式辅助分析

为补充分数之外的结构层面信息，本项目对每个 pocket 的复跑最优分子计算了 4 Å 内近邻残基。

| Pocket   | Top compound        |   复跑 score | 近邻残基摘要                                                                                               |
|:---------|:--------------------|-------------:|:-----------------------------------------------------------------------------------------------------------|
| BP1      | DB00776             |        -8.02 | LEU47B, LYS55B, ASP57B, SER50B, LEU56B, MET49B, ALA39B, ASN41B, LEU48B                                     |
| BP2      | 686393;CHEMBL392451 |       -10.76 | GLU185B, SER182B, ARG304B, MET169B, VAL186B, LEU188B, SER227B, ASN189B, ILE225B, PRO245B, GLY170B, ILE240B |

其中 BP1 命中残基包括 LEU47B、LYS55B、ASP57B、SER50B、MET49B、ALA39B、ASN41B，能对应 BP1 残基锚点；BP2 命中残基包括 GLU185B、ARG304B、MET169B、SER227B、ASN189B、ILE225B，也能对应 BP2 的关键区域。这个结果可用于结合模式讨论，但仍属于 docking pose 层面的计算证据。

![Binding-mode contact views](../results/figures/binding_mode_top_hits.png)

同时提供 `results/pymol/` 下的 PyMOL-ready 脚本；如果在本地恢复 receptor 和 docking pose 文件，可以用 PyMOL 进一步渲染 ray-traced 图或 2D interaction diagram。

## 图表输出

| 图 | 用途 |
|---|---|
| `results/figures/workflow_summary.png` | README 顶部流程与结果汇总图 |
| `results/figures/score_distribution_by_pocket.png` | BP1/BP2 原始分数分布对比 |
| `results/figures/top10_scores_by_pocket.png` | 各 pocket Top10 候选分子排序 |
| `results/figures/vina_rerun_vs_original_scores.png` | 原始分数与 Vina 复跑分数对比 |
| `results/figures/bp_score_scatter_overlap.png` | 两个 pocket 之间候选分子重叠与分数关系 |
| `results/figures/top_ligand_structures.png` | Top 候选分子 2D 结构图 |
| `results/figures/binding_mode_top_hits.png` | BP1/BP2 复跑最优分子的结合模式接触图 |
| `results/pymol/` | PyMOL-ready 结合模式渲染脚本 |
| `results/figures/enopt_style_reranking_vs_legacy.png` | EnOpt-style 加权分数与旧分数对比 |
| `results/figures/ensemble_score_matrix_heatmap.png` | compound × receptor conformation score matrix |
| `results/figures/enopt_style_conformation_weights.png` | 不同 receptor conformation 的权重 |

## EnOpt / ensemble 方法边界

参考文献中的 EnOpt 是用于提升 ensemble virtual screening 排序效果和可解释性的机器学习方法。严格意义上，EnOpt 通常需要：

- 多个蛋白构象或多个 docking score channels；
- 已知 active / decoy 或其他可监督标签；
- 可用于训练、验证和解释权重的 benchmark 数据。

本项目在原始 BP1/BP2 两口袋 score matrix 基础上，进一步生成了一个小型 FtsZ receptor ensemble：以重建后的 receptor 结构为起点，使用 ANM normal-mode perturbation 得到 base + mode1 ± + mode2 ± 共 5 个构象，并对 BP1/BP2 共 90 个 cleaned unique structures 进行多构象 Vina docking。当前 450 个 docking job 均完成，输出了 compound × conformation score matrix、ensemble-best / ensemble-mean 排序，以及 EnOpt-style conformation-weighted consensus score。

当前权重不是由实验活性标签监督训练得到，而是用旧 worksheet 排序的一致性作为探索性加权依据。因此对外建议表述为 **EnOpt-style weighted consensus reranking**，不要表述为严格 supervised EnOpt model。

关键输出包括：

- `results/tables/receptor_ensemble_manifest.csv`
- `results/tables/ensemble_vina_scores.csv`
- `results/tables/enopt_style_score_matrix.csv`
- `results/tables/enopt_style_weights.csv`
- `results/tables/enopt_style_top_hits.csv`
- `docs/enopt_extension.md`

当前 pilot ensemble 的加权 Top 结果：

| pocket | 加权 Top compound | weighted score, kcal/mol | legacy score, kcal/mol |
|---|---:|---:|---:|
| BP1 | CHEMBL329500 | -7.691 | -7.270 |
| BP2 | 85Z | -9.702 | -7.594 |

## 当前局限

- 原始 docking box/config 未找到，因此本项目 box 来自 pocket detection 输出和残基锚点的重建估计；
- 当前复跑覆盖清洗后的 90 个 unique structures；原始 Excel 的 method-level 记录已经合并、去盐和结构去重；
- 受体制备采用 OpenBabel PDBQT 流程；更严格的研究版本可进一步加入质子化状态、缺失残基、水分子和金属离子处理策略；
- docking score 仅是计算筛选指标，不等于 IC50、MIC 或其他实验活性；
- 当前 EnOpt-style 部分是 normal-mode receptor ensemble + weighted consensus pilot，不是最终 supervised ML model。

## 后续工作

1. 在现有结合模式接触图和 `results/pymol/` 脚本基础上，进一步生成 PyMOL ray-traced 图和 2D interaction diagram；
2. 用 curated experimental conformers 或 MD snapshots 替换当前 normal-mode pilot ensemble；
3. 寻找 FtsZ 已知 active / inactive 或 decoy set，训练并验证更严格的 supervised EnOpt-style reranking；
4. 将当前 pipeline 结构迁移到 GPCR 或其他靶点时，需要重新建立 receptor ensemble、小分子库抽样和验证设计。

## 可复现命令

公开仓库包含 processed tables 和最终 figures；如果要从原始 Excel / receptor / pocket files 完整复跑，需要先把原始输入文件按 `docs/data_availability.md` 放回 `data/raw/drive/`。

公开输出检查 / 图表再生成：

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m compileall -q scripts
python scripts/draw_top_ligands.py
python scripts/analyze_enopt_style.py
python scripts/make_summary_figure.py
```

本地完整复跑：

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
python scripts/draw_top_ligands.py
python scripts/draw_binding_mode_views.py
python scripts/make_summary_figure.py
python scripts/build_receptor_ensemble.py --target-ca-rmsd 0.65 --modes 2
python scripts/run_ensemble_vina.py --exhaustiveness 3 --cpu 4
python scripts/analyze_enopt_style.py
```
