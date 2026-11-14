"""空间—时间扫描的统计与搜索正确性测试（仅标准库）。"""
import copy
import itertools
import math
import random
import unittest
from collections import deque

from app.catalog import load
from app.scan import (
    ScanQuery, build_domain, find_clusters, llr, upper_bound_llr,
    simulate_null, iter_connected_sets, iter_windows, best_connected_set,
    _window_universe, inspect_scan, MODE_EXACT, MODE_PRUNED,
    DISJOINT_CELLS, DISJOINT_NODES, NullSimulation,
)
from app.service import SCHEMAS


def make_payload(nodes, edges, cells, data_version='t-1'):
    return {
        'schema_version': 1,
        'metadata': {'data_version': data_version, 'nodes': list(nodes)},
        'collections': {
            'cells': [
                {'id': f'c{i}', 'revision': 1, 'node': n, 'day': d,
                 'exposure': e, 'events': k}
                for i, (n, d, e, k) in enumerate(cells)
            ],
            'edges': [
                {'id': f'e{i}', 'revision': 1, 'from': a, 'to': b}
                for i, (a, b) in enumerate(edges)
            ],
        },
    }


def connected_brute_force(nodes, adjacency, subset):
    start = next(iter(subset))
    seen = {start}
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for w in adjacency[v]:
            if w in subset and w not in seen:
                seen.add(w)
                queue.append(w)
    return seen == set(subset)


# 6 节点路径图，3 天，节点 c 在中间三天制造持续超额
NODES = ['a', 'b', 'c', 'd', 'e', 'f']
EDGES = [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e'), ('e', 'f')]
CELLS = []
for d in range(3):
    for n in NODES:
        if n == 'c':
            CELLS.append((n, d, 100, 12 if d == 1 else 10))
        else:
            CELLS.append((n, d, 100, 2))
PAYLOAD = make_payload(NODES, EDGES, CELLS)


class EnumerationTests(unittest.TestCase):
    def test_connected_subsets_match_brute_force_without_duplicates(self):
        snap = load(PAYLOAD, SCHEMAS)
        domain = build_domain(snap)
        query = ScanQuery()
        lo, hi, universe = 0, 0, tuple(NODES)
        got = list(iter_connected_sets(domain, universe, set(), query))
        self.assertEqual(len(got), len(set(got)))
        adj = domain.adjacency
        expected = {
            frozenset(s) for r in range(1, len(NODES) + 1)
            for s in itertools.combinations(NODES, r)
            if connected_brute_force(NODES, adj, set(s))
        }
        self.assertEqual(set(got), expected)
        # 路径图连通子集数量为 n(n+1)/2 = 21
        self.assertEqual(len(got), 21)

    def test_size_cap_changes_candidate_domain_but_pruned_matches_exact(self):
        snap = load(PAYLOAD, SCHEMAS)
        domain = build_domain(snap)
        q = ScanQuery(max_spatial_size=2)
        lo, hi, universe = 0, 0, tuple(NODES)
        got = {s for s in iter_connected_sets(domain, universe, set(), q)}
        self.assertTrue(all(len(s) <= 2 for s in got))
        self.assertEqual(len(got), 11)  # 6 单点 + 5 条边
        exact = find_clusters(domain, ScanQuery(mode=MODE_EXACT, max_spatial_size=2))
        pruned = find_clusters(domain, ScanQuery(mode=MODE_PRUNED, max_spatial_size=2))
        self.assertEqual([c.candidate.key() for c in exact],
                         [c.candidate.key() for c in pruned])

    def test_missing_day_excludes_node_from_window_observable_domain(self):
        cells = [(n, d, 100, 1) for d in range(3) for n in NODES]
        cells = [c for c in cells if not (c[0] == 'c' and c[1] == 1)]  # c 缺一天
        snap = load(make_payload(NODES, EDGES, cells), SCHEMAS)
        domain = build_domain(snap)
        universe_wide = _window_universe(domain, 0, 2)
        self.assertNotIn('c', universe_wide)
        universe_narrow = _window_universe(domain, 0, 0)
        self.assertIn('c', universe_narrow)


class ScoreAndBoundTests(unittest.TestCase):
    def test_llr_basic_properties(self):
        # 非超额（事件份额 <= 暴露份额）截为 0
        self.assertEqual(llr(1, 10, 10, 100), 0.0)
        # 零事件、零暴露零事件均为 0，不产生无限值
        self.assertEqual(llr(0, 10, 0, 100), 0.0)
        self.assertEqual(llr(0, 0, 5, 100), 0.0)
        # 超额为正，且随事件数上升
        s1 = llr(3, 10, 10, 100)
        s2 = llr(4, 10, 10, 100)
        self.assertGreater(s1, 0.0)
        self.assertGreater(s2, s1)
        # 零暴露有事件是矛盾
        with self.assertRaisesRegex(ValueError, '矛盾'):
            llr(1, 0, 10, 100)
        # 与未截断 KL 形式在超额域一致
        n, mu, N, E = 4, 10, 10, 100
        p = mu / E
        kl = n * math.log(n / (N * p)) + (N - n) * math.log((N - n) / (N * (1 - p)))
        self.assertAlmostEqual(llr(n, mu, N, E), kl)

    def test_upper_bound_dominates_all_supersets_randomized(self):
        rng = random.Random(7)
        N, E = 40, 1000.0
        for trial in range(300):
            n = rng.randrange(0, N + 1)
            mu = rng.uniform(1.0, E * 0.9)
            addable = rng.randrange(0, N - n + 1)
            extra_mu = rng.uniform(0.0, E - mu)
            bound = upper_bound_llr(n, mu, addable, N, E)
            # 任意超集：n' ∈ [n, n+addable]，μ' >= μ
            actual = llr(n + addable, mu + extra_mu, N, E)
            self.assertGreaterEqual(bound + 1e-12, actual,
                                    f'上界失效 trial={trial}')
            # 角点本身（暴露不增、事件加满）取到上界
            self.assertAlmostEqual(bound, llr(n + addable, mu, N, E))
        # μ=0 无法有限上界时返回 inf，安全不剪枝
        self.assertEqual(upper_bound_llr(0, 0, 5, N, E), math.inf)

    def test_pruned_matches_exhaustive_on_random_graphs(self):
        rng = random.Random(11)
        for trial in range(25):
            n_nodes = rng.randint(4, 8)
            nodes = [str(i) for i in range(n_nodes)]
            # 随机连通图：先生成随机生成树再补边
            edges = set()
            order = nodes[:]
            rng.shuffle(order)
            for i in range(1, n_nodes):
                parent = order[rng.randrange(i)]
                child = order[i]
                self.assertNotEqual(parent, child)
                edges.add(tuple(sorted((parent, child))))
            for a, b in itertools.combinations(nodes, 2):
                if rng.random() < 0.25:
                    edges.add((a, b))
            days = rng.randint(1, 4)
            cells = []
            for d in range(days):
                for n in nodes:
                    e = rng.randint(20, 120)
                    k = rng.choices([0, 1, 2, 3, 4, 8], weights=[30, 30, 20, 10, 8, 2])[0]
                    cells.append((n, d, e, k))
            snap = load(make_payload(nodes, list(edges), cells), SCHEMAS)
            domain = build_domain(snap)
            kwargs_variants = [
                {},
                {'max_spatial_size': rng.randint(1, n_nodes)},
                {'max_exposure_share': rng.uniform(0.2, 0.8)},
                {'max_window_days': rng.randint(1, days)},
                {'disjoint_rule': DISJOINT_NODES},
            ]
            for kw in kwargs_variants:
                exact = find_clusters(domain, ScanQuery(mode=MODE_EXACT, **kw))
                pruned = find_clusters(domain, ScanQuery(mode=MODE_PRUNED, **kw))
                self.assertEqual(
                    [c.candidate.key() for c in exact],
                    [c.candidate.key() for c in pruned],
                    f'模式不一致 trial={trial} kw={kw}',
                )

    def test_best_connected_set_respects_banned_nodes(self):
        # 两个分离热点 a 与 f，禁入 a 区域后最优应落在 f 且节点不相交
        cells = []
        for d in range(3):
            for n in NODES:
                hot = n in ('a', 'f')
                cells.append((n, d, 100, 12 if hot else 2))
        snap = load(make_payload(NODES, EDGES, cells), SCHEMAS)
        domain = build_domain(snap)
        query = ScanQuery()
        lo, hi, universe = 0, 2, tuple(NODES)
        full = best_connected_set(domain, universe, lo, hi, set(), query)
        self.assertIsNotNone(full)
        banned = set(full.nodes)
        second = best_connected_set(domain, universe, lo, hi, banned, query)
        self.assertIsNotNone(second)
        self.assertFalse(set(second.nodes) & banned)


class ContradictionAndEdgeTests(unittest.TestCase):
    def test_zero_exposure_with_events_rejected(self):
        cells = [(n, d, 0 if (n == 'a' and d == 0) else 100,
                  1 if (n == 'a' and d == 0) else 0)
                 for d in range(2) for n in NODES]
        with self.assertRaisesRegex(ValueError, '矛盾'):
            build_domain(load(make_payload(NODES, EDGES, cells), SCHEMAS))

    def test_closed_days_zero_exposure_zero_events_do_not_inflate_risk(self):
        cells = []
        for d in range(3):
            for n in NODES:
                if n == 'a':
                    cells.append((n, d, 0, 0))  # 持续关闭
                else:
                    cells.append((n, d, 100, 1))
        snap = load(make_payload(NODES, EDGES, cells), SCHEMAS)
        domain = build_domain(snap)
        self.assertEqual(llr(0, 0, 5, domain.total_exposure), 0.0)
        result = inspect_scan(make_payload(NODES, EDGES, cells),
                              simulations=50, seed=3)
        for cluster in result['clusters']:
            self.assertNotIn('a', cluster['nodes'])
        self.assertTrue(math.isfinite(result['global_max_llr']))

    def test_no_events_yields_no_clusters_but_valid_report(self):
        cells = [(n, d, 100, 0) for d in range(2) for n in NODES]
        result = inspect_scan(make_payload(NODES, EDGES, cells), simulations=37)
        self.assertEqual(result['clusters'], [])
        self.assertEqual(result['global_max_llr'], 0.0)
        self.assertIsNone(result['global_p_value'])


class SimulationAndRevisionTests(unittest.TestCase):
    def test_p_values_and_max_correction_structure(self):
        result = inspect_scan(PAYLOAD, simulations=199, seed=5)
        self.assertEqual(len(result['null_simulation']['maxima']), 199)
        self.assertTrue(result['clusters'])
        top = result['clusters'][0]
        self.assertIn('c', top['nodes'])
        self.assertGreater(top['events'], top['expected_events'])
        self.assertGreater(top['llr'], 0.0)
        self.assertTrue(1 / 200 <= top['mc_p_value'] <= 1.0)
        # 全部发现共用同一次模拟的最大值分布
        for c in result['clusters']:
            self.assertTrue(c['mc_p_value'] >= result['global_p_value'] - 1e-12)

    def test_same_seed_reproduces_maxima(self):
        domain = build_domain(load(PAYLOAD, SCHEMAS))
        q = ScanQuery(simulations=80, seed=42)
        first = simulate_null(domain, q)
        second = simulate_null(domain, q)
        self.assertEqual(first.maxima, second.maxima)
        # 不同种子允许得到不同轨迹
        other = simulate_null(domain, ScanQuery(simulations=80, seed=43))
        self.assertNotEqual(first.maxima, other.maxima)

    def test_revision_invalidates_old_maximum_distribution(self):
        domain = build_domain(load(PAYLOAD, SCHEMAS))
        q = ScanQuery(simulations=60, seed=9)
        null_sim = simulate_null(domain, q)
        # 查询指纹变化（改上限）拒绝复用
        changed_query = ScanQuery(simulations=60, seed=9, max_spatial_size=3)
        with self.assertRaisesRegex(ValueError, '不可复用'):
            null_sim.p_value(1.0, domain.digest, changed_query.fingerprint())
        # 事件更正：摘要 digest 变化，旧分布拒绝复用
        revised = copy.deepcopy(PAYLOAD)
        row = next(r for r in revised['collections']['cells']
                   if r['node'] == 'c' and r['day'] == 1)
        row['revision'] = 2
        row['events'] = 20
        revised_domain = build_domain(load(revised, SCHEMAS))
        self.assertNotEqual(revised_domain.digest, domain.digest)
        with self.assertRaisesRegex(ValueError, '不可复用'):
            null_sim.p_value(1.0, revised_domain.digest, q.fingerprint())
        # 重新模拟后全局显著性可复算，且热点分数确实改变
        new_result = inspect_scan(revised, simulations=60, seed=9)
        old_result = inspect_scan(PAYLOAD, simulations=60, seed=9)
        self.assertNotEqual(new_result['global_max_llr'], old_result['global_max_llr'])
        self.assertEqual(new_result['input_digest'], revised_domain.digest)

    def test_closed_day_revision_changes_digest_and_recomputes(self):
        # 新增一条关闭日修订（零暴露零事件），digest 必须变化
        revised = copy.deepcopy(PAYLOAD)
        revised['collections']['cells'].append(
            {'id': 'a3', 'revision': 1, 'node': 'a', 'day': 3,
             'exposure': 0, 'events': 0})
        before = inspect_scan(PAYLOAD, simulations=50, seed=1)
        after = inspect_scan(revised, simulations=50, seed=1)
        self.assertNotEqual(before['input_digest'], after['input_digest'])
        self.assertEqual(after['observable_day_range'], [0, 3])
        # 同 seed 下模拟轨迹结构完整、可复现
        again = inspect_scan(revised, simulations=50, seed=1)
        self.assertEqual(after['null_simulation']['maxima'],
                         again['null_simulation']['maxima'])


class DisjointRuleTests(unittest.TestCase):
    def test_cell_disjoint_allows_same_node_at_nonoverlapping_times(self):
        cells = []
        for d in range(4):
            for n in NODES:
                burst = (n == 'a' and d in (0, 3))
                cells.append((n, d, 100, 15 if burst else 1))
        payload = make_payload(NODES, EDGES, cells)
        cell_rule = inspect_scan(payload, simulations=80, seed=2,
                                 disjoint_rule=DISJOINT_CELLS, max_findings=5)
        node_rule = inspect_scan(payload, simulations=80, seed=2,
                                 disjoint_rule=DISJOINT_NODES, max_findings=5)
        # cell 规则下允许同一节点出现在不相交时间段的两个发现中
        node_sets = [set(c['nodes']) for c in cell_rule['clusters']]
        windows = [(c['day_start'], c['day_end']) for c in cell_rule['clusters']]
        for (s1, w1), (s2, w2) in itertools.combinations(zip(node_sets, windows), 2):
            if s1 & s2:
                self.assertFalse(w1[0] <= w2[1] and w2[0] <= w1[1])
        # node 规则下所有发现节点集合两两不相交
        for s1, s2 in itertools.combinations(
                [set(c['nodes']) for c in node_rule['clusters']], 2):
            self.assertFalse(s1 & s2)
        # node 规则是更强约束，发现数不会更多
        self.assertLessEqual(len(node_rule['clusters']),
                             len(cell_rule['clusters']))


class ProvenanceAndContractTests(unittest.TestCase):
    def test_finding_trace_and_observable_domain_reported(self):
        result = inspect_scan(PAYLOAD, simulations=50, seed=6,
                              max_spatial_size=3)
        top = result['clusters'][0]
        # 可追溯到每个单元的事件与暴露、基线日期
        trace_keys = {(r['node'], r['day']) for r in top['cell_trace']}
        expected = {(n, d) for d in range(top['day_start'], top['day_end'] + 1)
                    for n in top['nodes']}
        self.assertEqual(trace_keys, expected)
        self.assertEqual(top['baseline_days'],
                         list(range(top['day_start'], top['day_end'] + 1)))
        self.assertEqual(top['window_observable_nodes'], NODES)
        trace_events = sum(r['events'] for r in top['cell_trace'])
        trace_exposure = sum(r['exposure'] for r in top['cell_trace'])
        self.assertEqual(trace_events, top['events'])
        self.assertEqual(trace_exposure, top['exposure'])
        # 候选域公开声明：非圆形、非仅局部最优
        cd = result['candidate_domain']
        self.assertFalse(cd['circular_scan'])
        self.assertFalse(cd['local_optima_only'])
        self.assertEqual(cd['max_spatial_size'], 3)

    def test_simulation_drawn_under_multinomial_null(self):
        # 零假设成立的合成数据：均匀暴露、事件率一致，多数运行不应极端显著
        cells = [(n, d, 100, 1) for d in range(4) for n in NODES]
        significant = 0
        for seed in range(12):
            result = inspect_scan(make_payload(NODES, EDGES, cells),
                                  simulations=199, seed=seed)
            if result['global_p_value'] is not None and result['global_p_value'] < 0.05:
                significant += 1
        self.assertLessEqual(significant, 2)  # 12 次中 0 或 1 次最可能，留宽限


class CandidateDomainTests(unittest.TestCase):
    def test_default_cap_is_disclosed_and_none_restores_full_domain(self):
        result = inspect_scan(PAYLOAD, simulations=30, seed=1)
        self.assertEqual(result['candidate_domain']['max_spatial_size'], 8)
        full = inspect_scan(PAYLOAD, simulations=30, seed=1, max_spatial_size=None)
        self.assertIsNone(full['candidate_domain']['max_spatial_size'])
        # 完整域不改变已选最优热点
        self.assertEqual(
            result['clusters'][0]['nodes'], full['clusters'][0]['nodes'])

    def test_exposure_share_cap_enforced_in_both_modes(self):
        for mode in (MODE_EXACT, MODE_PRUNED):
            result = inspect_scan(PAYLOAD, simulations=30, seed=1,
                                  mode=mode, max_exposure_share=0.3)
            for cluster in result['clusters']:
                self.assertLessEqual(
                    cluster['exposure'] / result['total_exposure'], 0.3 + 1e-12)

    def test_branch_and_bound_never_underestimates_optimum_randomized(self):
        # 直接在随机数据上对比每窗口 B&B 最优与完整穷举最优（无上限域）
        rng = random.Random(23)
        for trial in range(15):
            n_nodes = rng.randint(4, 7)
            nodes = [str(i) for i in range(n_nodes)]
            edges = set()
            order = nodes[:]; rng.shuffle(order)
            for i in range(1, n_nodes):
                edges.add(tuple(sorted((order[rng.randrange(i)], order[i]))))
            for a, b in itertools.combinations(nodes, 2):
                if rng.random() < 0.3:
                    edges.add((a, b))
            days = rng.randint(1, 3)
            cells = [(n, d, rng.randint(20, 120),
                      rng.choices([0, 1, 2, 3, 5], weights=[30, 30, 20, 12, 8])[0])
                     for d in range(days) for n in nodes]
            domain = build_domain(
                load(make_payload(nodes, list(edges), cells), SCHEMAS))
            q = ScanQuery(max_spatial_size=None)
            for lo, hi, universe in iter_windows(domain, q):
                got = best_connected_set(domain, universe, lo, hi, set(), q)
                ref = None
                for members in iter_connected_sets(domain, universe, set(), q):
                    n = sum(domain.node_events_in_window(v, lo, hi) for v in members)
                    mu = sum(domain.node_exposure_in_window(v, lo, hi) for v in members)
                    score = llr(n, mu, domain.total_events, domain.total_exposure)
                    if score <= 0:
                        continue
                    key = (-score, -n, lo, hi, tuple(sorted(members)))
                    if ref is None or key < ref:
                        ref = key
                if ref is None:
                    self.assertIsNone(got, f'trial={trial} window=({lo},{hi})')
                else:
                    self.assertEqual(got.key(), ref,
                                     f'剪枝低估最优 trial={trial} window=({lo},{hi})')


if __name__ == '__main__':
    unittest.main()
