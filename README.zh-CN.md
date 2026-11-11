# methylation-clock-eval

[English](README.md) · [简体中文](README.zh-CN.md)

> One-line summary: Cross-cohort evaluation of published DNA-methylation aging clocks (via `pyaging`) on two public whole-blood Illumina 450K cohorts (GSE40279, GSE87571), with explicit training-overlap flags, a resolved sample-count audit (732 GEO samples = 421 + 311, no repeat individuals), blood-cell-composition analysis of age acceleration, and an optional frozen ElasticNet external test.

这是一个**业余、纯计算（无湿实验）**的入门项目：在两个公开的全血 450K 甲基化队列上，复现并跨队列评估已发表的表观遗传时钟，重点是**不自欺**——记录训练重叠、缺失 CpG、偏差，而不是追求漂亮的数字。

> ⚠️ **非医疗建议，不是临床工具。** 只使用 GEO 上公开、已去标识的数据；这里**不做任何临床结论**，不能用于个人或临床的生物学年龄判断。时钟误差 ≠ 寿命证据，时钟"回拨" ≠ 真正变年轻。

## 数据

| 队列 | 来源 | 样本数 | 年龄 | 文件 |
|---|---|---|---|---|
| GSE40279 | Hannum et al. 2013, Mol Cell | 656（全部有年龄、性别） | 19–101 | series matrix（含 beta 值 + 元数据，470,043 个 cg 探针） |
| GSE87571 | Johansson et al. 2013, PLoS One（NSPHS 队列，瑞典北部；数据保存者单位为 Uppsala 大学） | 732（3 个缺年龄，2 个缺性别；评估用 729） | 14–94 | series matrix（仅元数据）+ 两个 supplementary beta 矩阵（每个样本的 `Xn` 列是 beta、`Xn.1` 列是检测 p 值，已核对） |

**732 vs 421 已解释**（见下文「GSE87571 样本数核查」）：GEO 的 732 条 = Johansson 2013 的 421 例 + 后续研究追加的 311 例，**没有发现同一个体的重复样本**；其中 3 条记录没有年龄/基因型/吸烟信息，不是可识别的个体，所以评估用 729 个个体。

## 目录结构

```
README.md                 本文件
requirements.txt          直接依赖（已安装的精确版本）；requirements.lock.txt 为完整 pip freeze
run_all.sh / Makefile     一键重跑
scripts/
  config.py               时钟列表、训练重叠标注表（含证据来源）、年龄分箱
  01_download.sh          下载 GEO 处理后矩阵（可续传，gzip -t 校验，写 sha256）
  02_clock_features.py    读取各时钟的元数据与 CpG 列表（权重来自 pyaging 的 Hugging Face 仓库）
  03_parse.py             流式解析：只保留时钟所需 CpG（约 2.2 万个）+ 元数据 + 全探针 QC
  04_run_clocks.py        用 pyaging 跑时钟；计算每个时钟/队列缺失 CpG 比例
  05_evaluate.py          指标、年龄分箱/性别分层残差、图
  06_train_elasticnet.py  （选做）在 GSE40279 训练 ElasticNetCV，冻结后外部测试 GSE87571
  07_sample_identity.py   GSE87571 的 732 vs 421：Sample_Group、基因型/甲基化重复个体检查、去重标记
  02b_cell_reference.py   用 Reinius 2012 纯化白细胞（GSE35069）构建 450K 细胞类型参考（仿 minfi pickCompProbes）
  08a_estimate_cells.py   Houseman 式约束投影估计 6 类白细胞比例（cell_deconv.py 为投影实现）
  08_cell_composition.py  年龄加速与细胞比例的相关、细胞校正后的比较与排名
  09_subset_metrics.py    GSE87571 去重/按 Sample_Group 分批的指标
data/            (gitignored) raw/ 与 processed/
results/         metrics.csv 等；plots/ 为 PNG
docs/PLAN.zh-CN.md     分析计划、已知坑、下一步
```

## 一键重跑

```bash
cd methylation-clock-eval
./run_all.sh                    # 已有 data/raw 时可用 SKIP_DOWNLOAD=1 ./run_all.sh ；或 make all
```
下载约 4 GB（含 Reinius 参考 GSE35069 约 120 MB）。在开发机（8 核 / 15 GB 内存、磁盘缓存已热）上，下载之后的全部分析步骤实测约 3 分钟；`pytest` 只跑轻量测试，不需要下载数据。Python 3.13，见 `requirements.txt`。

## 主要结果（均由本仓库代码实际算出；完整见 `results/metrics.csv`）

指标含义：MAE/RMSE 单位为年；bias = 预测 − 实际的均值（正 = 预测偏老）；`MAE_after_removing_bias` 仅为诊断（用评估标签去掉平均偏移后的 MAE）。
**训练重叠标记**：`yes` = 文献明确该队列在训练数据里（此时只是"复现检查"，不是独立验证）；`likely_yes`/`no_evidence`/`unknown` 的含义见 `scripts/config.py`。`no_evidence` 只表示没有找到训练包含的证据，**不等于证明独立**。

| 队列 | 时钟 | n | MAE | RMSE | R² | r | bias | 缺失CpG | 训练重叠 |
|---|---|---|---|---|---|---|---|---|---|
| GSE40279 | horvath2013 | 656 | 4.77 | 6.33 | 0.815 | 0.918 | -2.34 | 0/353 | yes |
| GSE40279 | hannum | 656 | 5.51 | 6.78 | 0.788 | 0.946 | +4.75 | 0/71 | yes |
| GSE40279 | skinandblood | 656 | 4.36 | 5.79 | 0.846 | 0.940 | -2.45 | 1/391 | unknown |
| GSE40279 | zhangen | 656 | 4.00 | 5.38 | 0.867 | 0.964 | -0.04 | 1/514 | likely_yes |
| GSE40279 | lin | 656 | 4.12 | 5.34 | 0.868 | 0.932 | -0.00 | 0/99 | yes |
| GSE40279 | intrinclock | 656 | 3.69 | 4.93 | 0.888 | 0.944 | -0.56 | 0/380 | unknown |
| GSE40279 | dnamphenoage* | 656 | 9.04 | 11.05 | 0.436 | 0.852 | -7.89 | 0/513 | no_evidence |
| GSE40279 | dunedinpace† | 656 | – | – | – | 0.135 | – | 16/20000 | no_evidence |
| GSE87571 | horvath2013 | 729 | 3.90 | 5.05 | 0.942 | 0.971 | +0.33 | 0/353 | no_evidence |
| GSE87571 | hannum | 729 | 5.13 | 6.14 | 0.914 | 0.975 | +4.01 | 0/71 | no_evidence |
| GSE87571 | skinandblood | 729 | 2.90 | 3.63 | 0.970 | 0.986 | +0.61 | 0/391 | unknown |
| GSE87571 | zhangen | 729 | 6.29 | 8.24 | 0.845 | 0.991 | -3.23 | 0/514 | no_evidence |
| GSE87571 | lin | 729 | 5.67 | 7.21 | 0.881 | 0.960 | +3.12 | 0/99 | no_evidence |
| GSE87571 | intrinclock | 729 | 3.21 | 4.17 | 0.960 | 0.982 | -1.32 | 0/380 | unknown |
| GSE87571 | dnamphenoage* | 729 | 9.34 | 10.78 | 0.735 | 0.968 | -8.90 | 0/513 | no_evidence |
| GSE87571 | dunedinpace† | 729 | – | – | – | 0.593 | – | 0/20000 | no_evidence |

\* DNAmPhenoAge 预测的是"表型年龄"（用死亡率相关临床指标训练），不是实足年龄，MAE/偏差在年龄尺度上没有意义，仅看相关性。
† DunedinPACE 预测的是衰老**速度**（EPIC 平台训练），不是年龄，所以不报 MAE：与实足年龄的 Pearson r 在 GSE40279 为 0.135、GSE87571 为 0.593。

选做：ElasticNetCV（`results/own_model_metrics.csv`；GSE40279 80/20 划分在任何特征处理之前，Imputer+Scaler 只在训练集上拟合，冻结后直接测试 GSE87571，不重新拟合）：

| 评估集 | n | MAE | RMSE | R² | r | bias |
|---|---|---|---|---|---|---|
| GSE40279 训练集（样本内） | 524 | 1.13 | 1.45 | 0.990 | 0.996 | 0.00 |
| GSE40279 留出 20% | 132 | 3.39 | 5.02 | 0.885 | 0.942 | +0.51 |
| GSE87571 外部（全部） | 729 | 2.85 | 3.65 | 0.970 | 0.986 | +0.37 |
| GSE87571 外部（≥19 岁） | 645 | 2.82 | 3.63 | 0.963 | 0.982 | +0.75 |

特征面板为"各时钟用到的 CpG 的并集中、两个队列都有的 21,739 个"（事先固定，没有用任何测试/外部标签挑选）；非零系数 581 个。注意该面板本身是发表研究挑选的，其中一些研究用过 GSE40279，所以留出集成绩偏乐观；外部成绩也受 GSE87571 年龄分布更宽（含 14–30 岁，年轻人多）影响，**不能和 GSE40279 上的成绩直接比较**。

## 关键观察与诚实的局限

1. **Hannum 在它自己的训练数据上也有 +4.7 年的系统偏差**（MAE 5.5，扣除偏差后 3.8）。pyaging 里 Hannum 的截距为 0（我检查了权重文件），所以输出是"原始系数直出"，没有校准；我没有去核对原论文的校准做法，也没有核对 GEO 矩阵与论文所用的归一化是否一致。所以这个复现**不能**说是"精确复现了论文的 3.9 年误差"。
2. 各时钟的偏差随年龄呈系统趋势（`results/residuals_by_age_bin.csv`）：Horvath/SkinAndBlood 对高龄者偏年轻（80+ 岁约 −3 到 −8 年），ZhangEN 在 30 岁以下严重偏年轻（约 −12 到 −14 年）。这是已知的"回归向均值/年龄范围"效应，不能简单解读为"生物学上更年轻"。
3. 排名会随队列变：GSE40279 上 IntrinClock 最好（MAE 3.7），GSE87571 上 SkinAndBlood 最好（2.9）；ZhangEN 的相关性最高（0.991）但 GSE87571 上 MAE 6.3，说明它的问题主要是校准（偏差和斜率），不是排序能力。
4. 数据处理上没有做 BMIQ/批次校正，直接用 GEO 已处理 beta 值（细胞组成只做了后验分析，见上文）；两个队列的预处理管线不同，一部分跨队列差异可能就来自这里。GSE87571 有少量 NA beta（每样本中位 0.12%，最多 3.3%）按 CpG 用队列均值填补；所有时钟所需 CpG 缺失比例都 ≤ 0.26%（见 `results/missing_cpgs.csv`），所以缺失 CpG 在这里不是主要误差来源。
5. 性别：Hannum 对男性偏老更多（GSE40279 男 +5.98 / 女 +3.59），其他时钟差异较小，详见 `results/residuals_by_sex.csv`。


## GSE87571 样本数核查（732 vs 421）

结论：**已解释，且没有重复个体。**

- GEO 的 `GSE87571_additional_sample_chararcteristics.xlsx`（保存者追加的样本特征）里有 `Sample_Group` 列：2012 组 **421** 个、2013 组 **311** 个，合计 732。2012 组的年龄中位数是 44，与 Johansson 2013 论文里的 N = 421、中位年龄 44、范围 14–94 吻合；所以 421 就是原论文那批，311 是后来追加的（`Sample_Group` 的含义是我根据数字吻合推断的，GEO 没有明说这一列的含义）。
- 追加的样本出自同一个队列（NSPHS，瑞典北部 Karesuando 2006 年与 Soppero 2009 年采样；`year_of_collection`：2006 年 683 个、2009 年 49 个）：Ahsan et al. 2017, PLoS Genetics 13(9):e1007005 用这套数据（GEO 登记号即 GSE87571），写明经质控后有 **729** 个样本有甲基化数据，另有 2 个遗传学离群样本被剔除、8 个对照样本（含 2 对重复）被移出。
- 我数据层面的核对：729 个有年龄且有基因型的个体；另有 **3 条记录没有年龄、基因型、吸烟信息**（GSM2334015、GSM2334342、GSM2334619），它们不是可识别的个体，已在评估中排除（这就是之前 729 vs 732 的差别）。它们是不是论文说的对照样本，**我没有证实**，只是与"对照样本"的描述一致。
- 重复个体检查：用 104 个基因型位点的剂量值做两两比较（267,546 对）：最小平均绝对差 **0.109**（中位 0.472，同一个体的重复应接近 0）；用 2.2 万个时钟相关 CpG 的中心化 beta 做相关：最大 **0.736**（中位 −0.003；技术重复预期接近 1）。设定阈值（基因型差 < 0.05 或甲基化相关 > 0.90）后找到的重复对为 **0**。基因型最接近的几对性别/年龄都不同，更像亲缘个体（这是人口较封闭的队列），不是重复样本；亲缘个体我没有去除，所以样本仍可能不完全独立。
- 去重评估：因为没有重复，去重后和之前的 729 例完全一致（见 `results/metrics_gse87571_subsets.csv`，`dedup_individuals` 与 `all_with_age` 逐项相同）。
- 论文提到的"同一个体相隔三年采血的重复样本"在这 729 个样本里没有出现（说明没有被放进 GEO 或已被剔除，我无法区分这两种可能）。

### 按批次分开看（`results/metrics_gse87571_subsets.csv`）

同一队列里，2012 组（421）与 2013 组（311）的偏差明显不同，例如 Hannum 偏差 +4.96 / +2.73，IntrinClock −0.28 / −2.72，ZhangEN −3.99 / −2.21，Lin +2.24 / +4.30（年）。两组的年龄分布也不同（中位 44 vs 53.5），而且没有校正批次，所以这里的差异**不能归因于单一原因**；但它说明时钟的"整体偏移"对批次/预处理很敏感。

## 血细胞组成分析

### 方法与局限

- **估计细胞比例**：Houseman 2012 约束投影（w ≥ 0，Σw ≤ 1，与 minfi 的 `projectCellType` 相同约束）。参考矩阵由 Reinius 2012 的 6 位男性供者纯化白细胞（GEO GSE35069，CD4T、CD8T、NK、B、单核、粒细胞）构建，选探针方式仿 minfi 的 `pickCompProbes`（F 检验 p < 1e-8，每类 50 个高甲基化 + 50 个低甲基化探针，并集共 **550** 个；GSE40279 中缺 10 个，GSE87571 无缺失）。
- **与保存者估计的对照**（GSE87571，保存者给出 CD8T/CD4T/NK/Bcell/Mono/Gran 估计，应该也是 Houseman 法，但我没有核实其参考与归一化）：我的估计与之的 Pearson r：Gran 0.99、NK 0.96、Bcell 0.93、Mono 0.92、CD4T 0.92、**CD8T 仅 0.66**（我的均值 0.040 对保存者 0.083，偏低一半）。所以 CD8T 估计质量较差，含 CD8T 的结论要谨慎。
- **未重新归一化**：minfi 会把参考与目标样本一起做分位数归一化，这里直接用 GEO 已发布的 beta 值（参考和两个队列的预处理管线都不同）。参考样本只有 6 人、全是男性，不是群体代表；我在 Reinius 的全血和 PBMC 样本（不在参考里）上做过检查：全血粒细胞估计 0.50，**PBMC 的粒细胞估计 0.12（应接近 0）**，说明投影有偏差；所以估计值适合做"相对差异/协变量"，不适合当成绝对的细胞比例。
- 没有用 pyaging 自带的 `DeconvoluteBloodEPIC*`：它的 600 个参考 CpG 里有 343–346 个（约 57%）在 450K 芯片上不存在，缺失值被默认值替代；它与保存者估计的相关明显更差（NK 0.69、CD8T 0.77、Mono 0.82，粒细胞均值 0.30 对 0.56），只作为诊断保留在 `results/cell_estimates_agreement_gse87571.csv`，不用于结论。
- 没有做更细的亚型（如初始/记忆 T 细胞），所以与 CD8T/CD4T 亚型相关的年龄变化仍在残差里。

### 定义（每队列、每个时钟）

- **AA（年龄加速）** = 对 `预测年龄 ~ a + b×实足年龄` 做 OLS 后的残差（去掉随年龄变化的系统偏差/压缩，所以与年龄线性正交）。
- **IEAA（内在式加速）** = 对 `AA ~ 5 类细胞比例` 回归后的残差（去掉第 6 类粒细胞，因为比例之和约为 1）。
- 校正细胞后，"相对实足年龄的 MAE"不再是合适目标；所以我比较的是：**RMSE 原始**（未校准）、**RMSE 重校准**（10 折交叉验证的 `年龄 ~ a + b×预测`，只去掉整体偏移和斜率）、**RMSE 重校准+细胞**（再加入 5 类细胞比例，看细胞组成能否改善还原年龄）、以及 **SD(AA)** 对 **SD(IEAA)**（校正细胞前后加速的离散度）。重校准用了评估队列自己的标签（交叉验证内部），是诊断性指标，不是可部署时钟。

### (a) 年龄加速与估计细胞比例的 Pearson 相关（`results/accel_cell_correlations.csv`）

GSE40279：

| 时钟 | CD8T | CD4T | NK | Bcell | Mono | Gran |
|---|---|---|---|---|---|---|
| horvath2013 | +0.11 | -0.18 | +0.07 | -0.08 | +0.09 | +0.05 |
| hannum | +0.04 | -0.33 | +0.05 | -0.15 | +0.08 | +0.21 |
| skinandblood | +0.01 | -0.07 | +0.07 | -0.09 | +0.06 | +0.02 |
| zhangen | -0.07 | -0.08 | -0.02 | -0.10 | +0.07 | +0.10 |
| lin | +0.10 | -0.18 | -0.02 | -0.09 | +0.07 | +0.10 |
| intrinclock | -0.07 | +0.05 | -0.05 | -0.00 | +0.00 | +0.02 |
| dnamphenoage | -0.11 | -0.27 | -0.08 | -0.14 | +0.15 | +0.27 |
| dunedinpace | -0.02 | -0.24 | -0.06 | -0.04 | +0.11 | +0.19 |

GSE87571（n=729）：

| 时钟 | CD8T | CD4T | NK | Bcell | Mono | Gran |
|---|---|---|---|---|---|---|
| horvath2013 | +0.18 | -0.13 | +0.24 | -0.08 | -0.16 | -0.09 |
| hannum | -0.04 | -0.29 | +0.23 | -0.16 | -0.01 | +0.12 |
| skinandblood | +0.11 | -0.01 | +0.18 | -0.12 | -0.13 | -0.09 |
| zhangen | +0.05 | -0.00 | +0.14 | -0.09 | -0.07 | -0.07 |
| lin | +0.13 | -0.05 | +0.15 | -0.06 | -0.13 | -0.07 |
| intrinclock | +0.02 | +0.01 | +0.01 | -0.01 | -0.11 | +0.01 |
| dnamphenoage | -0.01 | -0.26 | +0.02 | -0.18 | -0.06 | +0.24 |
| dunedinpace | -0.01 | -0.30 | -0.06 | -0.06 | +0.08 | +0.25 |

读法：
- **Hannum** 的加速与 CD4T 比例（−0.33 / −0.29）、粒细胞比例（+0.21 / +0.12）有稳定相关，两个队列方向一致；**DNAmPhenoAge**（CD4T −0.27/−0.26，粒细胞 +0.27/+0.24）和 **DunedinPACE**（CD4T −0.24/−0.30，粒细胞 +0.19/+0.25）也是如此。这与"这几个指标部分在读白细胞组成"的说法一致，但相关不等于因果，也受估计误差影响。
- **IntrinClock** 的相关最小（两个队列都 |r| ≤ 0.11），符合它"减少对免疫细胞组成敏感"的设计目标。
- Horvath、SkinAndBlood、Lin 在 GSE87571 里与 NK 呈 +0.15 到 +0.24 的相关，在 GSE40279 里只有约 0.07 或更低，**跨队列不一致**，不能当成稳定结论（GSE87571 里 NK 与年龄本身的相关也高达 0.42，而 GSE40279 里是 0.08，我没有查明原因，可能是队列差异或估计伪影）。

### (b) 细胞校正后的比较（均为 Reinius 参考估计；`results/accel_adjustment_summary.csv`）

| 队列 | 时钟 | SD(AA) | SD(IEAA) | 细胞解释AA方差 R² | RMSE原始 | RMSE重校准(CV) | RMSE重校准+细胞(CV) |
|---|---|---|---|---|---|---|---|
| GSE40279 | horvath2013 | 5.06 | 4.94 | 0.043 | 6.33 | 5.86 | 5.86 |
| GSE40279 | hannum | 4.24 | 3.97 | 0.125 | 6.78 | 4.78 | 4.63 |
| GSE40279 | skinandblood | 4.19 | 4.16 | 0.018 | 5.79 | 5.03 | 5.06 |
| GSE40279 | zhangen | 4.74 | 4.69 | 0.022 | 5.38 | 3.92 | 3.94 |
| GSE40279 | lin | 4.98 | 4.88 | 0.038 | 5.34 | 5.36 | 5.33 |
| GSE40279 | intrinclock | 4.33 | 4.32 | 0.007 | 4.93 | 4.86 | 4.83 |
| GSE40279 | dnamphenoage | 6.96 | 6.56 | 0.111 | 11.05 | 7.74 | 7.66 |
| GSE40279 | dunedinpace | 0.13 | 0.12 | 0.073 | - | - | - |
| GSE87571 | horvath2013 | 5.02 | 4.74 | 0.109 | 5.05 | 4.99 | 4.86 |
| GSE87571 | hannum | 4.50 | 4.18 | 0.140 | 6.14 | 4.66 | 4.47 |
| GSE87571 | skinandblood | 3.58 | 3.46 | 0.064 | 3.63 | 3.54 | 3.48 |
| GSE87571 | zhangen | 3.78 | 3.72 | 0.031 | 8.24 | 2.86 | 2.85 |
| GSE87571 | lin | 6.40 | 6.24 | 0.049 | 7.21 | 5.85 | 5.80 |
| GSE87571 | intrinclock | 3.93 | 3.90 | 0.011 | 4.17 | 3.96 | 3.95 |
| GSE87571 | dnamphenoage | 5.84 | 5.57 | 0.089 | 10.78 | 5.25 | 5.13 |
| GSE87571 | dunedinpace | 0.11 | 0.10 | 0.110 | - | - | - |

细胞比例对加速（AA）方差的解释比例（R²）：Hannum 12.5% / 14.0%，DNAmPhenoAge 11.1% / 8.9%，DunedinPACE 7.3% / 11.0%，Horvath 4.3% / 10.9%，其余多数 ≤ 6%，IntrinClock 最低（0.7% / 1.1%）。

排名（1 = 最好；仅 6 个实足年龄时钟；`results/clock_rankings_before_after_cells.csv`）：

| 队列 | 时钟 | 排名 RMSE原始 | RMSE重校准 | 重校准+细胞 | SD(AA) | SD(IEAA) |
|---|---|---|---|---|---|---|
| GSE40279 | horvath2013 | 5 | 6 | 6 | 6 | 6 |
| GSE40279 | hannum | 6 | 2 | 2 | 2 | 1 |
| GSE40279 | skinandblood | 4 | 4 | 4 | 1 | 2 |
| GSE40279 | zhangen | 3 | 1 | 1 | 4 | 4 |
| GSE40279 | lin | 2 | 5 | 5 | 5 | 5 |
| GSE40279 | intrinclock | 1 | 3 | 3 | 3 | 3 |
| GSE87571 | horvath2013 | 3 | 5 | 5 | 5 | 5 |
| GSE87571 | hannum | 4 | 4 | 4 | 4 | 4 |
| GSE87571 | skinandblood | 1 | 2 | 2 | 1 | 1 |
| GSE87571 | zhangen | 6 | 1 | 1 | 2 | 2 |
| GSE87571 | lin | 5 | 6 | 6 | 6 | 6 |
| GSE87571 | intrinclock | 2 | 3 | 3 | 3 | 3 |

- **排名变化主要来自"校准"，而不是细胞组成。** 原始 RMSE 排名与重校准后排名的 Kendall τ 只有 0.07（GSE40279）和 0.20（GSE87571）：例如 ZhangEN 在 GSE87571 上原始 RMSE 垫底（第 6），重校准后变成第 1（8.24 → 2.86 年）。
- **在重校准基础上再加细胞比例，排名完全不变**（τ = 1.00）；RMSE 改善很小（Hannum 4.78→4.63、4.66→4.47，其余多数不到 0.15 年）。
- 看加速的离散度：校正细胞前后排名只有一处互换（GSE40279 里 Hannum 与 SkinAndBlood 交换第 1/2 名），GSE87571 不变（SD 排名前后 Kendall τ：GSE40279 为 0.87，GSE87571 为 1.00）。
- 用保存者的 GSE87571 细胞估计做敏感性分析，结论一致（`cell_source = depositor_houseman` 行）。
- **适用范围**：以上只对"这两个队列、这套细胞估计"成立；细胞估计本身有偏差（见局限），所以"校正了细胞"不等于"校正了真实的细胞组成"。

## 输出文件

（`results/predictions.csv`、`cell_props_*.csv` 为逐样本表，本地重跑生成，已被 `.gitignore` 排除，不随仓库发布。）


- `results/metrics.csv` 每个时钟 × 队列：n、MAE、中位 AE、RMSE、R²、Pearson r、Spearman、平均带符号误差、缺失 CpG 比例、训练重叠标记与依据
- `results/residuals_by_age_bin.csv`、`results/residuals_by_sex.csv`
- `results/missing_cpgs.csv`、`results/pyaging_log_*.txt`（pyaging 自己的缺失特征警告）
- `results/qc_summary.csv`、`results/clock_info.csv`、`results/predictions.csv`
- `results/own_model_metrics.csv`
- `results/plots/pred_vs_age.png`、`residual_vs_age.png`、`dunedinpace_vs_age.png`
- 样本数核查：`results/gse87571_identity_summary.csv`、`gse87571_closest_genotype_pairs.csv`、`gse87571_duplicate_candidates.csv`（空）、`metrics_gse87571_subsets.csv`
- 细胞组成：`results/cell_reference_450k_reinius.csv`（参考矩阵）、`cell_props_houseman450k.csv`、`cell_props_summary.csv`、`cell_reference_coverage.csv`、`cell_method_sanity.csv`、`cell_estimates_agreement_gse87571.csv`、`accel_cell_correlations.csv`、`accel_adjustment_summary.csv`、`clock_rankings_before_after_cells.csv`
- 细胞图：`results/plots/accel_cell_corr_heatmap.png`、`rmse_raw_vs_recalibrated_vs_cells.png`

## 未做 / 未覆盖

PC 系列时钟（需要约 7.8 万个 CpG 和主成分投影）、GrimAge（需要性别+年龄输入）、ZhangBLUP（约 32 万 CpG）没有跑；细胞亚型、批次校正没做；亲缘个体没有去除。详见 `docs/PLAN.zh-CN.md`。

## 数据来源与引用

数据由脚本运行时从 GEO 下载，不随仓库分发；使用数据请引用原始论文。完整列表（含每个时钟的原论文，来自 pyaging 元数据）见 [`docs/REFERENCES.md`](docs/REFERENCES.md)：
GSE40279 — Hannum 等 2013（Mol Cell，doi:10.1016/j.molcel.2012.10.016）；GSE87571 — Johansson 等 2013（PLoS ONE，doi:10.1371/journal.pone.0067378）与 Ahsan 等 2017（PLoS Genet，doi:10.1371/journal.pgen.1007005）；
GSE35069 — Reinius 等 2012（PLoS ONE，doi:10.1371/journal.pone.0041361）；pyaging — de Lima Camillo 2024（Bioinformatics，doi:10.1093/bioinformatics/btae200）。
GEO 数据与 Reinius 参考保留各自的使用条款；MIT 许可证只覆盖本仓库的代码和文字。

## 声明

**非医疗建议，不是临床工具。** 时钟加速不是寿命、健康或任何干预效果的证据。引用本仓库见 `CITATION.cff`（作者名是占位）。
