"""现有业务摘要，不推断未观察到的业务关系。"""
def analyze(snapshot):
    seen, daily = set(), {}
    nodes = set(snapshot.metadata['nodes'])
    for row in snapshot.rows('cells'):
        key = row['node'], row['day']
        if key in seen or row['node'] not in nodes:
            raise ValueError('重复单元或未知节点')
        seen.add(key)
        if row['exposure'] < 0 or not isinstance(row['events'], int) or row['events'] < 0:
            raise ValueError('暴露或事件数无效')
        if row['exposure'] == 0 and row['events']:
            raise ValueError('零暴露单元有事件')
        group = daily.setdefault(row['day'], {'exposure': 0, 'events': 0, 'observed_nodes': []})
        group['exposure'] += row['exposure']
        group['events'] += row['events']
        group['observed_nodes'].append(row['node'])
    return {str(k): {**v, 'observed_nodes': sorted(v['observed_nodes'])} for k, v in sorted(daily.items())}
