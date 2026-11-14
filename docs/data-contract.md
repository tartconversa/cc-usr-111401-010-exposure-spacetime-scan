# 资料与现有计算约定

## 修订快照

schema_version为整数1，metadata为业务方案字典，collections为业务集合字典。每条记录有非空字符串id和正整数revision。
同集合相同id取最高revision；同id同revision内容不同即冲突，即便存在更高修订也拒绝。deleted=true为撤回，之后更高修订可恢复。
撤回记录不要求携带业务字段。输入顺序与相同记录重试不改变语义；非有限数字拒绝。
Snapshot.rows返回独立副本；input_digest绑定有效修订和全部metadata。输出data_version沿用metadata，不由程序自造资料版本。
既有摘要入口为app.service.inspect(payload)，输出schema_version、input_digest、data_version和summary。

## 业务资料

cells 的node为图节点、day为整数日号、exposure为非负暴露量、events为非负整数；同node/day只有一个有效单元，未提供单元表示资料缺失，不等同于关闭。edges是无向邻接；nodes明确空间全集。baseline.analyze 按日汇总完整单元的暴露及事件，不生成热点。只有具备候选所需全部日单元的区域参与该窗口扫描；查询结果须记录此可观察域。

各集合必要字段：{"cells": ["node", "day", "exposure", "events"], "edges": ["from", "to"]}。

## 当前能力

app.baseline.analyze(snapshot)产生已有摘要，app.service.inspect负责绑定输入身份；新增分析可以增设独立入口，已有接口及上述数据语义保持兼容。

## 空间—时间扫描（app.scan）

入口 app.scan.inspect_scan(payload, **查询参数) / app.scan.run(snapshot, **查询参数)，CLI：`python -m app.scan 资料JSON [模拟次数]`。

- 候选域：邻接图上的全部连通子集 × 全部连续整数日段。节点在窗口内缺任一日子单元即不可观察，逐窗口从可观察域排除并记录到 window_observable_nodes；缺测不等同关闭。查询参数 max_spatial_size、max_exposure_share（候选暴露占总暴露 E 的份额上限）、max_window_days 收窄候选域，全部在输出 candidate_domain 公开；默认 max_spatial_size=8 是为大图可终止的受限搜索，小图传 None 恢复完整域。本模块不做圆形扫描、不只保留局部最优（candidate_domain 显式声明 circular_scan/local_optima_only 均为 false）。
- 零假设与分数：查询内基线暴露固定，条件于总事件数 N，事件按单元暴露比例 p_c=e_c/E 多项分配；候选归组为二项，分数 llr() 为似然比且仅超额取正，与 Poisson 扫描条件化后代数等价。零暴露零事件记 0 分（关闭日不产生无限风险）；零暴露却有事件在 build_domain 即拒绝（资料矛盾）。
- 搜索：mode=exact 枚举全部连通子集（小图穷举参照）；mode=pruned 用可证上界的分支限界，上界 LLR(n+n_R, μ) 中 n_R 仅统计前沿在 U\\excluded 子图内可达节点（有规模上限时再按预算取事件最多的若干节点），严格小于当前最优才剪枝，候选域与 exact 相同。随机图测试断言两模式结果一致。
- 最大值校正：simulate_null 每次模拟在同一候选域重扫并只保留全体候选最大 LLR；p_value=(1+不小于观测的模拟数)/(R+1)，主热点与次级发现共用同一最大分布。NullSimulation 绑定 input_digest 与 query_fingerprint，任一不匹配即拒绝给结论；事件更正或关闭日修订改变 digest，必须重新模拟，不得复用旧最大值分布。
- 次级发现：disjoint_rule=cells（默认，不共享重叠时段的节点单元，允许同节点异时段）或 nodes（任意时间节点不相交）。每个发现输出 nodes、起止日、baseline_days、events、exposure、expected_events、relative_risk、llr、mc_p_value、逐日 daily_baseline 与逐单元 cell_trace，可回溯事件、暴露与基线日期。
