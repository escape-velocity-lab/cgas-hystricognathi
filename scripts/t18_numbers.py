# -*- coding: utf-8 -*-
"""任务 #18b：改写稿件需要的新数字，全部现取并存档（坑 43）。

输出 results/t18_numbers.json 与 t18_numbers.txt。
"""
import json
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

NMR = {"431": "S", "479": "E", "495": "Y", "498": "T"}


def main():
    d = pd.read_csv(os.path.join(RES, "t13_species_sites.csv"))
    chen = pd.read_csv(os.path.join(RES, "t18_chen_figS3K.csv"))

    n = {}
    n["chen_n_species"] = len(chen)
    n["chen_overlap_with_us"] = int((chen["本项目_431"] != "未收录").sum())
    n["chen_sites_agree"] = int((chen["结论"] == "一致").sum())
    n["chen_4of4"] = chen.loc[chen["Chen_裸鼹鼠型命中数"] == 4, "图注物种名"].tolist()
    n["chen_3of4"] = chen.loc[chen["Chen_裸鼹鼠型命中数"] == 3, "图注物种名"].tolist()
    n["chen_hystricognath_rows"] = chen.loc[chen["Chen_分类群标注"] == "Ctenohystrica", "图注物种名"].tolist()
    n["chen_has_cavia"] = bool(chen["拉丁名"].str.startswith("Cavia").any())
    n["chen_lifespan_vs_anage"] = {
        r["图注物种名"]: [r["Chen_寿命标注_年"], r["本项目_最大寿命_年"]]
        for _, r in chen.iterrows()
        if str(r["本项目_最大寿命_年"]) != "" and abs(float(r["本项目_最大寿命_年"]) - float(r["Chen_寿命标注_年"])) > 1
    }

    # 啮齿目外的分布（对应 Chen 图注 "None of these changes are observed in long-lived
    # species outside Rodentia"）
    nonrod = d[d["目"] != "Rodentia"]
    n["nonrodent_taxa"] = int(len(nonrod))
    for s, aa in NMR.items():
        sub = nonrod[nonrod[s] == aa]
        n["nonrodent_site%s_nmr_n" % s] = int(len(sub))
        n["nonrodent_site%s_nmr_species" % s] = sub["物种"].tolist()
    n["nonrodent_max_hits"] = int(nonrod["命中"].max())
    n["nonrodent_hits_ge2"] = int((nonrod["命中"] >= 2).sum())
    # 啮齿目外的长寿物种（LQ >= 2）里有没有命中
    ll = nonrod[(nonrod["LQ"] >= 2)]
    n["nonrodent_LQge2_n"] = int(len(ll))
    n["nonrodent_LQge2_with_hit"] = ll.loc[ll["命中"] >= 1, "物种"].tolist()

    # 豪猪亚目（含豚鼠）
    hy = d[d["科"].isin([
        "Bathyergidae", "Heterocephalidae", "Caviidae", "Chinchillidae", "Octodontidae",
        "Ctenomyidae", "Echimyidae", "Dasyproctidae", "Cuniculidae", "Hystricidae",
        "Erethizontidae", "Dinomyidae", "Petromuridae", "Thryonomyidae", "Myocastoridae", "Abrocomidae",
    ])]
    n["hyst_taxa"] = int(len(hy))
    n["hyst_4of4"] = hy.loc[hy["命中"] == 4, "物种"].tolist()
    n["hyst_3of4"] = hy.loc[hy["命中"] == 3, "物种"].tolist()
    n["hyst_families"] = sorted(hy["科"].unique().tolist())

    # Chen 的 18 种里本项目有 LQ 的，LQ 排名
    n["cavia_rows"] = d.loc[d["物种"].str.startswith("Cavia"),
                            ["物种", "431", "479", "495", "498", "最大寿命_年", "LQ", "命中"]].to_dict("records")
    n["degu_chinchilla"] = d.loc[d["物种"].isin(["Octodon degus", "Chinchilla lanigera"]),
                                 ["物种", "最大寿命_年", "LQ", "命中"]].to_dict("records")

    with open(os.path.join(RES, "t18_numbers.json"), "w", encoding="utf-8") as f:
        json.dump(n, f, ensure_ascii=False, indent=2, default=str)
    with open(os.path.join(RES, "t18_numbers.txt"), "w", encoding="utf-8") as f:
        for k, v in n.items():
            f.write("%s: %s\n" % (k, v))
    for k, v in n.items():
        print("%s: %s" % (k, v))


if __name__ == "__main__":
    main()
