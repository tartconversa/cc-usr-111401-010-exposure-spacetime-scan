# 乡村活动拥挤热点的暴露校正与联合扫描

本项目提供空间连通与连续时间候选相关资料的现有业务摘要。当前代码执行资料修订读取与基础统计，不包含扩展分析。

app/catalog.py 按集合和记录键解析修订，拒绝同修订冲突；app/baseline.py 承载当前业务计算；app/service.py 输出摘要及输入摘要哈希。
app/scan.py 提供暴露校正的空间—时间扫描：连通子集×连续日段候选、条件多项零假设的似然比分、带可证上界的分支限界（小图另有穷举参照）、同次模拟最大值校正的全局显著性，以及绑定资料/查询指纹的零分布复算保护。统计约定见[资料契约](docs/data-contract.md)。
数据字段及现有接口见 [资料契约](docs/data-contract.md)。fixtures/sample.json 是合成小数据，不对应真实村庄或个人。
运行环境为 Python 3.14，现有模块仅使用标准库。

读取并计算当前摘要：

    python -m app.service fixtures/sample.json

运行空间—时间扫描（可选第二参数为蒙特卡洛模拟次数）：

    python -m app.scan fixtures/sample.json 999

运行当前模块检查：

    python -m unittest discover -s tests -v

检查源码编译：

    python -m compileall -q app tests
