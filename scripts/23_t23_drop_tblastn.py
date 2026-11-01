"""23（任务 #23d-b，来自 #21 自审「需补」第 1 条）: 剔除只有 tblastn 片段的物种后重算全树 FEL。
动机：117 个密码子（人 406-522）恰好是低质量 tblastn 序列唯一贡献数据的区段。比对误差会抬高 dN、
把 ω 推向 1 —— 正好是我们阴性结论的方向。把这些物种全部去掉再跑一次，看 4 个位点的结论是否不变。
剔除标准：results/t17_codon/sources.csv 里「来源」以 efetch 开头（= 任务 #13 的 tblastn 重建），不是 RefSeq/GenBank CDS。
输出：results/t17_codon/aln_cds_only.fasta、tree_cds_only.nwk
用法: python 23_t23_drop_tblastn.py <proj_root>
"""
import sys, os, io
import pandas as pd
from Bio import SeqIO, Phylo

ROOT = sys.argv[1]
R = lambda *p: os.path.join(ROOT, *p)
src = pd.read_csv(R('results', 't17_codon', 'sources.csv'), encoding='utf-8-sig')
src.columns = ['species', 'order', 'family', 'source', 'cov', 'hystri']
drop = set(src.loc[src['source'].str.startswith('efetch'), 'species'].str.replace(' ', '_'))
print('剔除 %d 个 tblastn 重建物种：' % len(drop))
for s in sorted(drop):
    print('   ', s)

recs = list(SeqIO.parse(R('results', 't17_codon', 'aln.fasta'), 'fasta'))
keep = [r for r in recs if r.id not in drop]
print('比对：%d -> %d 条' % (len(recs), len(keep)))
assert len(recs) - len(keep) == len(drop), '名字对不上：%s' % (drop - {r.id for r in recs})
out_aln = R('results', 't17_codon', 'aln_cds_only.fasta')
SeqIO.write(keep, out_aln, 'fasta')

tree = Phylo.read(R('results', 't17_codon', 'tree_plain.nwk'), 'newick')
for name in drop:
    tree.prune(name)
# 修掉剪枝后可能出现的单子节点
tree.collapse_all(lambda c: c.is_preterminal() is False and len(c.clades) == 1)
tips = {t.name for t in tree.get_terminals()}
print('树：剪枝后 %d 个叶子' % len(tips))
assert tips == {r.id for r in keep}, '树和比对的物种集不一致'
out_tree = R('results', 't17_codon', 'tree_cds_only.nwk')
Phylo.write(tree, out_tree, 'newick')
print('写出:', out_aln); print('写出:', out_tree)
