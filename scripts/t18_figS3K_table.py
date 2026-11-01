# -*- coding: utf-8 -*-
"""任务 #18a：把 Chen et al. 2025 Science 补充材料 Fig. S3K 的 18 物种比对整理成表，
并与本项目（#11/#13）自己算出的 4 位点残基逐格核对。

Fig. S3K 的两段比对在裸鼹鼠坐标下分别是 461-471 和 510-533；
标注的 4 个位点 463 / 511 / 527 / 530 对应人 cGAS 编号 431 / 479 / 495 / 498。
比对行文本从 PDF 第 5 页提取（PyMuPDF），并与 300 dpi 渲染图逐行目视核对。
"""
import csv
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

# (图注物种名, 图注寿命标注, 拉丁名, block1 461-471, block2 510-533, 图上的分类群标注)
ROWS = [
    ("Naked mole-rat", 37, "Heterocephalus glaber", "ELSKFCSYHVK", "EELEHYFIPGVNLFSQEYIDTKSK", "Ctenohystrica"),
    ("Damaraland mole-rat", 20, "Fukomys damarensis", "ELSEFCSYHMK", "EELEHYFIPGVNFFSQEYIDTKCK", "Ctenohystrica"),
    ("Chinchilla", 17, "Chinchilla lanigera", "ELSKFCSYHVK", "EELKHFFIPEVNLFSPEYINKRSK", "Ctenohystrica"),
    ("Degu", 8, "Octodon degus", "EMSKFCSYHVK", "EELKHFFIPEVNLFSQEYIDQKSK", "Ctenohystrica"),
    ("Alpine marmot", 17, "Marmota marmota", "ELDKFCSYHVK", "EQLEHYFIPVVNLFSRNHIDKISK", "Sciuromorpha"),
    ("Woodchuck", 14, "Marmota monax", "ELDKFCSYHVK", "EQLEHYFIPVVNLFSRNHIDKISK", "Sciuromorpha"),
    ("Gray squirrel", 23, "Sciurus carolinensis", "ELDKFCSYHVK", "EQLEHYFVPIVNLFSRKYIDKISK", "Sciuromorpha"),
    ("Rat", 5, "Rattus norvegicus", "ELDAFCSYHVK", "EKLDHYFIPKFNLFSQELIDKKSK", "Mouse-related clade"),
    ("Mouse", 4, "Mus musculus", "ELDAFCSYHVK", "EKLDHYFIPKFNLFSQELIDRKSK", "Mouse-related clade"),
    ("Blind mole-rat", 21, "Nannospalax galili", "ELDKFCSYHVK", "EKLEHFFIPKFNLFSEELVDTKSK", "Mouse-related clade"),
    ("Ord's kangaroo rat", 9, "Dipodomys ordii", "ELNKFCSYHVK", "ENLAHYFIPVVNLFSQDKIDNVSR", "Mouse-related clade"),
    ("Greater horseshoe bat", 30, "Rhinolophus ferrumequinum", "HLDDFCSYHVK", "EDLKHYFIPGFNLFSRDLIDRASK", "Other mammalian species"),
    ("African elephant", 65, "Loxodonta africana", "ELNKFCSYHVK", "EQLRHYFIPGFNLFSPDEVNKISK", "Other mammalian species"),
    ("Asian elephant", 79, "Elephas maximus", "ELNKFCSYHVK", "EQLRHYFIPGFNLFSPDEVNKISK", "Other mammalian species"),
    ("Blue whale", 110, "Balaenoptera musculus", "ELDKFCSYHVK", "EQLEHYFIPEVNLFSRDDIDKPSK", "Other mammalian species"),
    ("Sperm whale", 77, "Physeter macrocephalus", "ELGKFCSYHVK", "EQLEHYFIPEVNLFSRDHIDKLSK", "Other mammalian species"),
    ("Rhesus monkey", 40, "Macaca mulatta", "HLDKFSSYHVK", "EQLKNYFIPGFNLFSSNLIDKTSK", "Other mammalian species"),
    ("Human", 122, "Homo sapiens", "HLDKFSSYHVK", "EKLENYFIPEFNLFSSNLIDKRSK", "Other mammalian species"),
]

NMR = {"463": "S", "511": "E", "527": "Y", "530": "T"}
# 裸鼹鼠坐标 → 人坐标
MAP = {"463": "431", "511": "479", "527": "495", "530": "498"}

B1_START = 461
B2_START = 510


def read_site(b1, b2, pos):
    p = int(pos)
    if p < B2_START:
        return b1[p - B1_START]
    return b2[p - B2_START]


def main():
    ours = pd.read_csv(os.path.join(RES, "t13_species_sites.csv"))
    ours = ours.set_index("物种")

    out = []
    for name, yrs, latin, b1, b2, clade in ROWS:
        rec = {
            "图注物种名": name,
            "拉丁名": latin,
            "Chen_寿命标注_年": yrs,
            "Chen_分类群标注": clade,
        }
        hit = 0
        agree = 0
        n_cmp = 0
        for pos in ("463", "511", "527", "530"):
            aa = read_site(b1, b2, pos)
            rec["Chen_%s" % pos] = aa
            if aa == NMR[pos]:
                hit += 1
            hp = MAP[pos]
            if latin in ours.index:
                mine = str(ours.loc[latin, hp])
                rec["本项目_%s" % hp] = mine
                n_cmp += 1
                if mine == aa:
                    agree += 1
            else:
                rec["本项目_%s" % hp] = "未收录"
        rec["Chen_裸鼹鼠型命中数"] = hit
        if latin in ours.index:
            rec["本项目_裸鼹鼠型命中数"] = int(ours.loc[latin, "命中"])
            rec["本项目_最大寿命_年"] = ours.loc[latin, "最大寿命_年"]
            rec["本项目_LQ"] = ours.loc[latin, "LQ"]
            rec["本项目_序列来源"] = ours.loc[latin, "来源"]
        else:
            rec["本项目_裸鼹鼠型命中数"] = ""
            rec["本项目_最大寿命_年"] = ""
            rec["本项目_LQ"] = ""
            rec["本项目_序列来源"] = ""
        rec["逐位点一致"] = "%d/%d" % (agree, n_cmp) if n_cmp else "无法比较"
        rec["结论"] = "一致" if n_cmp and agree == n_cmp else ("不一致" if n_cmp else "本项目未收录")
        out.append(rec)

    cols = ["图注物种名", "拉丁名", "Chen_寿命标注_年", "Chen_分类群标注",
            "Chen_463", "Chen_511", "Chen_527", "Chen_530", "Chen_裸鼹鼠型命中数",
            "本项目_431", "本项目_479", "本项目_495", "本项目_498", "本项目_裸鼹鼠型命中数",
            "本项目_最大寿命_年", "本项目_LQ", "本项目_序列来源", "逐位点一致", "结论"]
    df = pd.DataFrame(out)[cols]
    path = os.path.join(RES, "t18_chen_figS3K.csv")
    df.to_csv(path, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL)
    print(df.to_string())
    print("\n写出:", path)

    # 关键核查
    print("\n[核查] 豚鼠属是否出现在 Fig. S3K:", any("Cavia" in r[2] for r in ROWS))
    print("[核查] 全部 18 行逐位点一致:", (df["结论"] == "一致").sum(), "/", len(df))
    print("[核查] Chen 图中 4/4 的物种:", df.loc[df["Chen_裸鼹鼠型命中数"] == 4, "图注物种名"].tolist())
    print("[核查] Chen 图中 3/4 的物种:", df.loc[df["Chen_裸鼹鼠型命中数"] == 3, "图注物种名"].tolist())
    print("[核查] 寿命标注差异（Chen vs 本项目 AnAge，差 >1 年）:")
    for _, r in df.iterrows():
        if r["本项目_最大寿命_年"] != "" and abs(float(r["本项目_最大寿命_年"]) - r["Chen_寿命标注_年"]) > 1:
            print("   ", r["图注物种名"], r["Chen_寿命标注_年"], "vs", r["本项目_最大寿命_年"])


if __name__ == "__main__":
    main()
