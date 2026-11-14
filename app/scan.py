"""空间—时间扫描统计：连通子集 × 连续日段候选上的暴露校正热点检测。

候选域（不在此处偷偷收窄）：
  空间：项目邻接图上的全部连通子集，可限制最大节点数与候选暴露份额 μ/E；
  时间：观测日号范围内的全部连续整数日段，可限制最大长度；
  可观察域：某节点在窗口内缺任一日子单元即不可观察，不参与该窗口扫描，
  结果逐窗口记录实际可观察节点（缺测不等同于关闭）。

零假设：本次查询内暴露基线固定，条件于总事件数 N，事件按各单元暴露比例
p_c = e_c/E 服从多项分布；候选归组后为二项分布，扫描分数为二项似然比
（与 Poisson 扫描统计在条件化后代数等价），仅超额取正：
  LLR = n·ln(n/(Np)) + (N−n)·ln((N−n)/(N(1−p)))，p = μ/E。
零暴露零事件（如关闭日）记 0 分，不制造无限风险；零暴露却有事件属资料矛盾。

显著性：蒙特卡洛在同一候选域上重复模拟，每次只保留全体候选的最大分数，
用同一次模拟的最大值分布校正全部候选的搜索多重性，而不是逐候选普通检验。
NullSimulation 绑定 input_digest 与查询指纹；事件更正或关闭日修订后摘要
变化、digest 变化，旧最大值分布即不可复用，必须重新模拟。

两种搜索模式候选域完全相同：
  exact  —— 小图穷举参照，枚举全部连通子集后按分数贪心选取不相交发现；
  pruned —— 分支限界，上界 LLR(n+n_R, μ)（n_R 为分支内剩余可纳入节点的
            事件合计）可证不低于任一超集分数，仅在严格小于当前最优时剪枝，
            次级发现通过逐轮加入不相交约束后重新求最优得到，与穷举贪心等价。
本模块不做圆形扫描、不做只保留局部最优的收窄；如未来改为受限候选域，
须在输出 candidate_domain 中公开改变内容。
"""
from __future__ import annotations

import hashlib
import math
import random
from collections import deque
from dataclasses import dataclass

from .catalog import canonical
from .service import SCHEMAS, load

MODE_EXACT = 'exact'
MODE_PRUNED = 'pruned'
DISJOINT_CELLS = 'cells'      # 不共享任何 (节点, 日) 单元
DISJOINT_NODES = 'nodes'      # 任意时间都不共享节点


# ---------------------------------------------------------------- 基础校验

@dataclass(frozen=True)
class ScanQuery:
    # 大图连通子集呈指数增长，默认限制空间规模使搜索可终止；
    # 小图参照可显式传 None 恢复完整候选域（exact 穷举或无上限 B&B）。
    max_spatial_size: int | None = 8
    max_exposure_share: float | None = None
    max_window_days: int | None = None
    max_findings: int = 10
    simulations: int = 999
    seed: int | None = 20261114
    disjoint_rule: str = DISJOINT_CELLS
    mode: str = MODE_PRUNED

    def fingerprint(self):
        payload = {
            'max_spatial_size': self.max_spatial_size,
            'max_exposure_share': self.max_exposure_share,
            'max_window_days': self.max_window_days,
            'max_findings': self.max_findings,
            'simulations': self.simulations,
            'seed': self.seed,
            'disjoint_rule': self.disjoint_rule,
            'mode': self.mode,
        }
        return hashlib.sha256(canonical(payload).encode('utf-8')).hexdigest()

    def validate(self):
        if self.mode not in (MODE_EXACT, MODE_PRUNED):
            raise ValueError('未知搜索模式')
        if self.disjoint_rule not in (DISJOINT_CELLS, DISJOINT_NODES):
            raise ValueError('未知不相交规则')
        if self.max_spatial_size is not None and self.max_spatial_size < 1:
            raise ValueError('最大空间规模至少为1')
        if self.max_window_days is not None and self.max_window_days < 1:
            raise ValueError('最大窗口长度至少为1')
        if not isinstance(self.max_findings, int) or self.max_findings < 1:
            raise ValueError('最大发现数至少为1')
        if not isinstance(self.simulations, int) or self.simulations < 1:
            raise ValueError('模拟次数至少为1')
        if self.seed is not None and not isinstance(self.seed, int):
            raise ValueError('随机种子应为整数')
        if self.max_exposure_share is not None and not (
            0.0 < self.max_exposure_share <= 1.0
        ):
            raise ValueError('暴露份额上限应在 (0, 1] 内')


@dataclass(frozen=True)
class GridDomain:
    digest: str
    nodes: tuple[str, ...]
    adjacency: dict[str, frozenset[str]]
    days: tuple[int, ...]
    cells: dict[tuple[str, int], tuple[float, int]]
    total_exposure: float
    total_events: int

    def node_events_in_window(self, node, day_start, day_end, event_source=None):
        source = self.cells if event_source is None else event_source
        return sum(source[(node, d)][1] for d in range(day_start, day_end + 1))

    def node_exposure_in_window(self, node, day_start, day_end):
        return sum(self.cells[(node, d)][0] for d in range(day_start, day_end + 1))


def build_domain(snapshot):
    """校验修订快照并构造扫描网格；零暴露有事件等矛盾在此拒绝。"""
    nodes = tuple(sorted(snapshot.metadata['nodes']))
    node_set = set(nodes)
    if not nodes:
        raise ValueError('空间节点全集为空')
    adjacency = {n: set() for n in nodes}
    for row in snapshot.rows('edges'):
        a, b = row['from'], row['to']
        if a not in node_set or b not in node_set:
            raise ValueError('邻接边包含未知节点')
        if a == b:
            raise ValueError('不接受自环边')
        adjacency[a].add(b)
        adjacency[b].add(a)

    cells: dict[tuple[str, int], tuple[float, int]] = {}
    seen = set()
    total_exposure = 0.0
    total_events = 0
    for row in snapshot.rows('cells'):
        node, day = row['node'], row['day']
        if node not in node_set:
            raise ValueError('单元包含未知节点')
        if type(day) is not int or isinstance(day, bool):
            raise ValueError('日号应为整数')
        key = node, day
        if key in seen:
            raise ValueError('重复单元')
        seen.add(key)
        exposure, events = row['exposure'], row['events']
        if not isinstance(exposure, (int, float)) or isinstance(exposure, bool) or exposure < 0:
            raise ValueError('暴露量无效')
        if type(events) is not int or events < 0:
            raise ValueError('事件数应为非负整数')
        if exposure == 0 and events > 0:
            raise ValueError('资料矛盾：零暴露单元却有事件')
        cells[key] = float(exposure), events
        total_exposure += float(exposure)
        total_events += events
    if not cells:
        raise ValueError('无任何观测单元，无法扫描')
    if total_events > 0 and total_exposure <= 0:
        raise ValueError('资料矛盾：存在事件但总暴露为零')
    days = tuple(sorted({d for _, d in cells}))
    return GridDomain(
        digest=snapshot.digest,
        nodes=nodes,
        adjacency={n: frozenset(nb) for n, nb in adjacency.items()},
        days=days,
        cells=cells,
        total_exposure=total_exposure,
        total_events=total_events,
    )


# ---------------------------------------------------------------- 分数

def _xlogx(x):
    return 0.0 if x == 0 else x * math.log(x)


def llr(events, exposure, total_events, total_exposure):
    """条件多项零假设下的二项似然比；非超额返回 0，零暴露零事件返回 0。"""
    if exposure < 0 or events < 0:
        raise ValueError('计数不能为负')
    if exposure == 0:
        if events > 0:
            raise ValueError('资料矛盾：零暴露候选却有事件')
        return 0.0
    if total_events == 0 or events == 0:
        return 0.0
    if events * total_exposure <= total_events * exposure:
        return 0.0  # 事件份额不超过暴露份额，非超额
    p = exposure / total_exposure
    value = 0.0
    if events:
        value += events * math.log(events / (total_events * p))
    outside = total_events - events
    if outside:
        value += outside * math.log(outside / (total_events * (1.0 - p)))
    value = max(0.0, value)
    return 0.0 if value == 0.0 else value  # 归一化 -0.0，保证比较一致


def upper_bound_llr(events, exposure, addable_events, total_events, total_exposure):
    """分支内任一超集的 LLR 上界。

    固定 N、E 时 LLR 对 n 单调不减、对 μ 单调不增（零假设侧被截为 0 仍保持
    单调），故超集取“事件加满、暴露不增”的角点即支配所有可行 (n', μ')。
    μ=0 时角点落入奇异域，无法给出有限上界，返回 +inf（安全地不剪枝）。
    """
    if exposure == 0:
        return math.inf
    return llr(events + addable_events, exposure, total_events, total_exposure)


# ---------------------------------------------------------------- 候选枚举

def _window_universe(domain, day_start, day_end):
    """窗口内每天都有单元的节点；缺测节点排除出该窗口的可观察域。"""
    width = day_end - day_start + 1
    counts = {n: 0 for n in domain.nodes}
    for d in range(day_start, day_end + 1):
        for n in domain.nodes:
            if (n, d) in domain.cells:
                counts[n] += 1
    return tuple(n for n in domain.nodes if counts[n] == width)


def iter_windows(domain, query):
    dmin, dmax = domain.days[0], domain.days[-1]
    for lo in range(dmin, dmax + 1):
        hi_cap = dmax if query.max_window_days is None else min(dmax, lo + query.max_window_days - 1)
        for hi in range(lo, hi_cap + 1):
            universe = _window_universe(domain, lo, hi)
            if universe:
                yield lo, hi, universe


def _initial_excluded(universe, root_index, banned):
    return set(universe[:root_index]) | banned


def iter_connected_sets(domain, universe, banned, query):
    """穷举连通子集：每个集合仅在其最小序号节点为根时生成一次。

    规范 include/exclude DFS：每个前沿节点 v 要么纳入并把其邻居并入前沿，
    要么永久排除（v 的邻居若可经其他路径到达，仍会留在/进入前沿），
    两分支互斥穷尽，同一连通集合恰好生成一次。
    """
    us = set(universe)
    banned = banned or set()
    size_cap = query.max_spatial_size

    def grow(members, frontier, excluded):
        yield frozenset(members)
        if size_cap is not None and len(members) >= size_cap or not frontier:
            return
        pending = set(frontier)
        while pending:
            v = max(pending)
            pending.remove(v)
            added = {
                w for w in domain.adjacency[v]
                if w in us and w not in excluded and w not in members and w not in pending
            }
            # 副本：排除集合沿各自分支独立累积，不污染兄弟分支
            yield from grow(members | {v}, pending | added, set(excluded))
            excluded = excluded | {v}  # 永久排除：排除分支不再回到 v

    for i, root in enumerate(universe):
        if root in banned:
            continue
        excluded = _initial_excluded(universe, i, banned)
        frontier = {
            v for v in domain.adjacency[root]
            if v in us and v not in excluded
        }
        yield from grow(frozenset((root,)), frontier, excluded)


@dataclass
class _Candidate:
    nodes: tuple[str, ...]
    day_start: int
    day_end: int
    events: int
    exposure: float
    score: float

    def key(self):
        # 确定性并列排序：分数、事件数、时间位置、空间构成
        return (-self.score, -self.events, self.day_start, self.day_end, self.nodes)


def _window_node_totals(domain, universe, day_start, day_end, event_source=None):
    source = domain.cells if event_source is None else event_source
    ev, ex = {}, {}
    for n in universe:
        ev[n] = sum(source[(n, d)][1] for d in range(day_start, day_end + 1))
        ex[n] = sum(domain.cells[(n, d)][0] for d in range(day_start, day_end + 1))
    return ev, ex


def _share_ok(exposure, domain, query):
    if query.max_exposure_share is None:
        return True
    return exposure / domain.total_exposure <= query.max_exposure_share


def best_connected_set(domain, universe, day_start, day_end, banned, query,
                       event_source=None):
    """分支限界求窗口内最优连通候选；剪枝只删上界严格低于 incumbent 的分支。

    event_source 给定时（蒙特卡洛模拟），事件计数取自模拟网格，暴露仍用基线。
    """
    us = set(universe)
    banned = banned or set()
    node_events, node_exposure = _window_node_totals(
        domain, universe, day_start, day_end, event_source
    )
    size_cap = query.max_spatial_size
    N, E = domain.total_events, domain.total_exposure
    best: _Candidate | None = None

    def consider(members, n, mu):
        nonlocal best
        if not _share_ok(mu, domain, query):
            return  # 超集暴露只会更大，调用方负责剪枝
        score = llr(n, mu, N, E)
        if score <= 0.0:
            return
        cand = _Candidate(tuple(sorted(members)), day_start, day_end, n, mu, score)
        if best is None or cand.key() < best.key():
            best = cand

    def _reachable(frontier, members, excluded):
        """分支内任何连通超集可能纳入的节点：在 U\\excluded 子图内、
        由前沿出发的可达闭包（excluded 中的永久排除点不可穿过）。"""
        seen = set(frontier)
        queue = deque(frontier)
        while queue:
            v = queue.popleft()
            for w in domain.adjacency[v]:
                if (w in us and w not in excluded and w not in members
                        and w not in seen):
                    seen.add(w)
                    queue.append(w)
        return seen

    def _addable_events(frontier, members, excluded, root_events):
        reachable = _reachable(frontier, members, excluded)
        events = [node_events[v] for v in reachable]
        if size_cap is not None:
            # 超集至多再纳入 (size_cap - |members|) 个节点，取事件最多的若干个
            budget = size_cap - len(members)
            if budget <= 0:
                return root_events
            events = sorted(events, reverse=True)[:budget]
        return root_events + sum(events)

    def grow(members, frontier, excluded, n, mu):
        nonlocal best
        if members and _share_ok(mu, domain, query):
            consider(members, n, mu)
        if not _share_ok(mu, domain, query):
            return  # 份额已超限，所有超集同样超限
        if size_cap is not None and len(members) >= size_cap or not frontier:
            return
        addable = _addable_events(frontier, members, excluded, n)
        bound = upper_bound_llr(n, mu, addable - n, N, E)
        if best is not None and bound < best.score:
            return  # 严格小于：并列时保留分支，结论仍精确
        pending = set(frontier)
        while pending:
            v = max(pending)
            pending.remove(v)
            added = {
                w for w in domain.adjacency[v]
                if w in us and w not in excluded and w not in members and w not in pending
            }
            # 副本：排除集合沿各自分支独立累积，与穷举枚举同一候选遍历顺序
            grow(members | {v}, pending | added, set(excluded),
                 n + node_events[v], mu + node_exposure[v])
            excluded = excluded | {v}

    # 热启动：每个允许根做一次贪心扩张，得到域内合法 incumbent，
    # 让剪枝在穷举早期即生效；不改变候选域与精确性。
    def warm_start(root, initial_excluded):
        members = frozenset((root,))
        n, mu = node_events[root], node_exposure[root]
        consider(members, n, mu)
        frontier = {
            v for v in domain.adjacency[root]
            if v in us and v not in initial_excluded
        }
        while frontier:
            if size_cap is not None and len(members) >= size_cap:
                break
            best_v = max(frontier, key=lambda v: node_events[v])
            mu_next = mu + node_exposure[best_v]
            if not _share_ok(mu_next, domain, query):
                frontier.discard(best_v)
                continue
            members = members | {best_v}
            n, mu = n + node_events[best_v], mu_next
            consider(members, n, mu)
            frontier.discard(best_v)
            frontier |= {
                w for w in domain.adjacency[best_v]
                if w in us and w not in initial_excluded
                and w not in members and w not in frontier
            }

    for i, root in enumerate(universe):
        if root in banned:
            continue
        warm_start(root, _initial_excluded(universe, i, banned))
    for i, root in enumerate(universe):
        if root in banned:
            continue
        excluded = _initial_excluded(universe, i, banned)
        frontier = {
            v for v in domain.adjacency[root]
            if v in us and v not in excluded
        }
        grow(frozenset((root,)), frontier, excluded,
             node_events[root], node_exposure[root])
    return best


# ---------------------------------------------------------------- 不相交选择

def _intervals_overlap(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


def _conflicts(candidate, selected_nodes, selected_windows, rule):
    own_nodes = set(candidate.nodes)
    for nodes, (w_start, w_end) in zip(selected_nodes, selected_windows):
        if not (own_nodes & nodes):
            continue
        if rule == DISJOINT_NODES:
            return True
        if _intervals_overlap(candidate.day_start, candidate.day_end, w_start, w_end):
            return True
    return False


def _banned_for_window(selected_nodes, selected_windows, day_start, day_end, rule):
    banned = set()
    for nodes, (w_start, w_end) in zip(selected_nodes, selected_windows):
        if rule == DISJOINT_NODES or _intervals_overlap(day_start, day_end, w_start, w_end):
            banned |= nodes
    return banned


# ---------------------------------------------------------------- 观测扫描

@dataclass
class _Found:
    candidate: _Candidate
    observable_nodes: tuple[str, ...]


def _exhaustive_scan(domain, query):
    """小图参照：枚举全部候选，按确定性顺序贪心保留不相交发现。"""
    candidates = []
    observability = {}
    for lo, hi, universe in iter_windows(domain, query):
        observability[(lo, hi)] = universe
        for members in iter_connected_sets(domain, universe, set(), query):
            n = sum(domain.node_events_in_window(v, lo, hi) for v in members)
            mu = sum(domain.node_exposure_in_window(v, lo, hi) for v in members)
            if not _share_ok(mu, domain, query):
                continue
            score = llr(n, mu, domain.total_events, domain.total_exposure)
            if score > 0.0:
                candidates.append(_Candidate(tuple(sorted(members)), lo, hi, n, mu, score))
    candidates.sort(key=lambda c: c.key())
    selected_nodes, selected_windows, picked = [], [], []
    for cand in candidates:
        if _conflicts(cand, selected_nodes, selected_windows, query.disjoint_rule):
            continue
        picked.append(_Found(cand, observability[(cand.day_start, cand.day_end)]))
        selected_nodes.append(set(cand.nodes))
        selected_windows.append((cand.day_start, cand.day_end))
        if len(picked) >= query.max_findings:
            break
    return picked


def _pruned_scan(domain, query):
    """逐轮求最优：把与已选发现冲突的节点作为禁入约束传给分支限界。

    cell 规则下仅在时间重叠的窗口禁入对应节点（不重叠窗口同节点仍可纳入，
    等价于 _conflicts 的逐节点形式）；node 规则下所有窗口全局禁入。
    B&B 在受约束候选域上求精确最优，与穷举模式的贪心选择结果一致。
    """
    picked = []
    selected_nodes, selected_windows = [], []
    windows = list(iter_windows(domain, query))
    while len(picked) < query.max_findings:
        incumbent = None
        incumbent_universe = None
        for lo, hi, universe in windows:
            banned = _banned_for_window(
                selected_nodes, selected_windows, lo, hi, query.disjoint_rule
            )
            cand = best_connected_set(domain, universe, lo, hi, banned, query)
            if cand is not None and (incumbent is None or cand.key() < incumbent.key()):
                incumbent, incumbent_universe = cand, universe
        if incumbent is None:
            break
        picked.append(_Found(incumbent, incumbent_universe))
        selected_nodes.append(set(incumbent.nodes))
        selected_windows.append((incumbent.day_start, incumbent.day_end))
    return picked


def find_clusters(domain, query):
    query.validate()
    if query.mode == MODE_EXACT:
        return _exhaustive_scan(domain, query)
    return _pruned_scan(domain, query)


def maximum_score(domain, query, event_source=None):
    """一次（模拟）扫描中全体候选的最大分数；蒙特卡洛最大值校正用。"""
    best = 0.0
    for lo, hi, universe in iter_windows(domain, query):
        cand = best_connected_set(
            domain, universe, lo, hi, set(), query, event_source=event_source
        )
        if cand is not None and cand.score > best:
            best = cand.score
    return best


# ---------------------------------------------------------------- 蒙特卡洛

@dataclass(frozen=True)
class NullSimulation:
    input_digest: str
    query_fingerprint: str
    seed: int
    simulations: int
    maxima: tuple[float, ...]

    def p_value(self, score, digest, fingerprint):
        """仅当资料摘要与查询指纹均未变化时，该最大值分布才有效。"""
        if digest != self.input_digest or fingerprint != self.query_fingerprint:
            raise ValueError('资料或查询已修订，旧最大值分布不可复用，请重新模拟')
        greater_equal = 1 + sum(1 for t in self.maxima if t >= score)
        return greater_equal / (self.simulations + 1)


def _simulate_multinomial_events(domain, rng):
    """条件于总事件数，按暴露比例多项分配事件到各单元。"""
    keys = sorted(domain.cells)
    weights = [domain.cells[k][0] for k in keys]
    draws = rng.choices(keys, weights=weights, k=domain.total_events)
    counts = {k: 0 for k in keys}
    for k in draws:
        counts[k] += 1
    return {k: (domain.cells[k][0], counts[k]) for k in keys}


def simulate_null(domain, query):
    """每次调用都重新模拟；本函数不提供跨摘要复用的缓存。"""
    query.validate()
    if domain.total_events == 0:
        return NullSimulation(domain.digest, query.fingerprint(),
                              query.seed if query.seed is not None else 0,
                              query.simulations, ())
    rng = random.Random(query.seed)
    maxima = []
    for _ in range(query.simulations):
        simulated = _simulate_multinomial_events(domain, rng)
        maxima.append(maximum_score(domain, query, event_source=simulated))
    return NullSimulation(
        input_digest=domain.digest,
        query_fingerprint=query.fingerprint(),
        seed=query.seed if query.seed is not None else 0,
        simulations=query.simulations,
        maxima=tuple(maxima),
    )


# ---------------------------------------------------------------- 结果组装

def _trace_cells(domain, found):
    c = found.candidate
    rows = []
    for d in range(c.day_start, c.day_end + 1):
        for node in c.nodes:
            exposure, events = domain.cells[(node, d)]
            rows.append({'node': node, 'day': d,
                         'exposure': exposure, 'events': events})
    return rows


def _daily_baseline(domain, found):
    c = found.candidate
    rows = []
    for d in range(c.day_start, c.day_end + 1):
        exposure = sum(domain.cells[(n, d)][0] for n in c.nodes)
        events = sum(domain.cells[(n, d)][1] for n in c.nodes)
        rows.append({'day': d, 'exposure': exposure, 'events': events})
    return rows


def _candidate_domain_description(query):
    return {
        'spatial': 'all_connected_subsets_on_project_adjacency_graph',
        'temporal': 'all_contiguous_integer_day_intervals',
        'max_spatial_size': query.max_spatial_size,
        'max_exposure_share': query.max_exposure_share,
        'max_window_days': query.max_window_days,
        'observability': 'nodes_missing_any_day_cell_are_excluded_per_window',
        'search_mode': query.mode,
        'pruning': (
            'none_exhaustive_reference' if query.mode == MODE_EXACT
            else 'branch_and_bound_with_proven_upper_bound_same_candidate_domain'
        ),
        'circular_scan': False,
        'local_optima_only': False,
    }


def run(snapshot, **params):
    """完整扫描：构造域、找热点、重新模拟全局零分布并绑定资料身份。"""
    query = ScanQuery(**params)
    query.validate()
    domain = build_domain(snapshot)
    findings = find_clusters(domain, query)
    null_sim = simulate_null(domain, query)
    fingerprint = query.fingerprint()

    clusters = []
    for found in findings:
        c = found.candidate
        expected = domain.total_events * c.exposure / domain.total_exposure
        relative_risk = (
            (c.events / c.exposure) / (domain.total_events / domain.total_exposure)
            if c.exposure > 0 else None
        )
        p_value = (
            null_sim.p_value(c.score, domain.digest, fingerprint)
            if domain.total_events > 0 else None
        )
        clusters.append({
            'nodes': list(c.nodes),
            'day_start': c.day_start,
            'day_end': c.day_end,
            'baseline_days': list(range(c.day_start, c.day_end + 1)),
            'events': c.events,
            'exposure': c.exposure,
            'expected_events': expected,
            'relative_risk': relative_risk,
            'llr': c.score,
            'mc_p_value': p_value,
            'window_observable_nodes': list(found.observable_nodes),
            'daily_baseline': _daily_baseline(domain, found),
            'cell_trace': _trace_cells(domain, found),
        })

    observed_max = max((f.candidate.score for f in findings), default=0.0)
    return {
        'schema_version': 1,
        'input_digest': domain.digest,
        'data_version': snapshot.metadata.get('data_version'),
        'query_fingerprint': fingerprint,
        'candidate_domain': _candidate_domain_description(query),
        'disjoint_rule': query.disjoint_rule,
        'observable_day_range': [domain.days[0], domain.days[-1]],
        'total_events': domain.total_events,
        'total_exposure': domain.total_exposure,
        'global_max_llr': observed_max,
        'global_p_value': (
            null_sim.p_value(observed_max, domain.digest, fingerprint)
            if domain.total_events > 0 and findings else None
        ),
        'null_simulation': {
            'seed': null_sim.seed,
            'simulations': null_sim.simulations,
            'maxima': list(null_sim.maxima),
            'note': '仅对 input_digest 与 query_fingerprint 同时有效；修订后必须重算',
        },
        'contradictions': [],
        'clusters': clusters,
    }


def inspect_scan(payload, **params):
    """从原始修订资料出发的便捷入口。"""
    snapshot = load(payload, SCHEMAS)
    return run(snapshot, **params)


def main():
    import json
    import sys
    from pathlib import Path

    if len(sys.argv) not in (2, 3):
        raise SystemExit('用法：python -m app.scan 资料JSON [模拟次数]')
    payload = json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
    params = {}
    if len(sys.argv) == 3:
        params['simulations'] = int(sys.argv[2])
    print(json.dumps(inspect_scan(payload, **params),
                     ensure_ascii=False, allow_nan=False, indent=2))


if __name__ == '__main__':
    main()
