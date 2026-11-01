"""30（任务 #25c）：把 draft_v3.md 转成投稿用的 Word 稿。

做法：
  1. 读 draft_v3.md，删掉开头 1~3 行的内部工作注释（不能跟着稿子投出去）
  2. 把 Figure legends 一节拆开，图片 + 图注嵌到正文对应位置（Fig.1 → §4 末、Fig.2 → §6 末、Fig.S1 → 补充材料）
  3. 给三张没有题注的表格补上 Table 题注
  4. pandoc → docx
  5. python-docx 后处理：Times New Roman 12 pt、双倍行距、表格加框线、全文行号（投稿常要求）

用法: python src/30_t25_export_docx.py [--zh 中文 md] [--read 阅读版] [--out 输出 docx]
解释器: 用 projects/01-epigenetic-clock/.venv（python-docx 装在那里；02 的 .venv 没有，见坑 86）
  投稿版: ../01-epigenetic-clock/.venv/Scripts/python.exe src/30_t25_export_docx.py
  阅读版: ... src/30_t25_export_docx.py --read --out manuscript/draft_v3_阅读版.docx
"""
import os
import re
import subprocess
import sys

import docx
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt, Inches, Cm, RGBColor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
PANDOC = os.environ.get('PANDOC', 'pandoc')

ZH = '--zh' in sys.argv
READ = '--read' in sys.argv   # 任务 #33：阅读版（1.35 倍行距、无行号、两端对齐、深蓝分级标题、斑马表格）


def _opt_path(flag):
    """--zh / --src 后面跟的路径（没跟路径就返回 None）。任务 #36：
    docstring 一直写的是 `--zh 中文 md`，但原来只把 --zh 当开关，源文件仍走 --src，
    按文档写的命令会去读英文稿然后报「找不到图注一节」。两种写法现在都支持。"""
    if flag not in sys.argv:
        return None
    i = sys.argv.index(flag) + 1
    if i < len(sys.argv) and not sys.argv[i].startswith('--'):
        return sys.argv[i]
    return None


SRC = _opt_path('--src') or _opt_path('--zh') or R('manuscript', 'draft_v3.md')
OUT = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else R('manuscript', 'draft_v3.docx')
TMP = os.path.splitext(OUT)[0] + '_export.md'

FIGDIR = R('manuscript', 'figures')

# 图注所在节的标题（中英两版）
LEGEND_HEAD = '## 图注' if ZH else '## Figure legends'
SUPP_HEAD = '## 补充材料' if ZH else '## Supplementary material'


def build_export_md(src, dst):
    text = open(src, encoding='utf-8').read()
    # 1. 删内部注释
    text = re.sub(r'^<!--.*?-->\s*', '', text, flags=re.S)

    lines = text.split('\n')
    # 2. 取出图注
    legends = {}
    try:
        i0 = lines.index(LEGEND_HEAD)
    except ValueError:
        raise SystemExit(f'找不到图注一节: {LEGEND_HEAD}')
    i1 = lines.index(SUPP_HEAD)
    for ln in lines[i0 + 1:i1]:
        m = re.match(r'\*\*(Figure|图) ?(S?\d)[.．]', ln.strip())
        if m:
            legends[m.group(2)] = ln.strip()
    assert set(legends) == {'1', '2', 'S1'}, legends
    # 删掉整节（含上一行可能的分隔）
    del lines[i0:i1]

    def fig_block(key, fname):
        return ['', f'![]({FIGDIR}/{fname}){{width=6.2in}}', '', legends[key], '']

    out, used = [], set()
    for ln in lines:
        out.append(ln)
        # Fig.1 放到 §4（祖先重建）末尾，Fig.2 放到 §6（结构）末尾
        if ln.startswith('### ') and ('5.' in ln[:8] or '五、' in ln[:8]) and '1' not in used:
            out[-1:] = fig_block('1', 'Figure1.png') + [ln]
            used.add('1')
        if ln.startswith('### ') and ('7.' in ln[:8] or '七、' in ln[:8]) and '2' not in used:
            out[-1:] = fig_block('2', 'Figure2.png') + [ln]
            used.add('2')
        if ln == SUPP_HEAD and 'S1' not in used:
            out += fig_block('S1', 'FigureS1.png')
            used.add('S1')
    assert used == {'1', '2', 'S1'}, used

    text = '\n'.join(out)
    # 3. 表题注：两版 md 都已自带，按出现顺序编号（1=§3 寿命、2=§5 选择、3=§7 候选）。
    # 2026-10-22 任务 #35 改了英文分支，任务 #36 改了中文分支：这里再注入会重复，
    # 且曾产生 2->3->1 的错误顺序。
    pass
    open(dst, 'w', encoding='utf-8').write(text)
    return dst


def run_pandoc(md, out):
    cmd = [PANDOC, md, '-f', 'markdown+pipe_tables+tex_math_dollars', '-o', out,
           '--standalone', '--resource-path', FIGDIR]
    subprocess.run(cmd, check=True)


def polish(path):
    d = docx.Document(path)
    font_name = 'Times New Roman'
    cjk = '宋体' if ZH else None

    def set_font(style):
        style.font.name = font_name
        style.font.size = Pt(12)
        rpr = style.element.get_or_add_rPr()
        rf = rpr.find(qn('w:rFonts'))
        if rf is None:
            rf = rpr.makeelement(qn('w:rFonts'), {})
            rpr.append(rf)
        rf.set(qn('w:ascii'), font_name)
        rf.set(qn('w:hAnsi'), font_name)
        if cjk:
            rf.set(qn('w:eastAsia'), cjk)

    normal = d.styles['Normal']
    set_font(normal)
    normal.paragraph_format.line_spacing = 2.0
    normal.paragraph_format.space_after = Pt(6)

    for name, size in [('Title', 16), ('Heading 1', 14), ('Heading 2', 12), ('Heading 3', 12)]:
        try:
            st = d.styles[name]
        except KeyError:
            continue
        set_font(st)
        st.font.size = Pt(size)
        st.font.bold = True
        st.paragraph_format.line_spacing = 2.0

    # 表格：加框线、表内单倍行距、字号 10
    for t in d.tables:
        try:
            t.style = d.styles['Table Grid']
        except KeyError:
            pass
        for row in t.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    p.paragraph_format.line_spacing = 1.0
                    p.paragraph_format.space_after = Pt(0)
                    for r in p.runs:
                        r.font.size = Pt(10)
                        r.font.name = font_name

    # 图片居中、限宽
    for p in d.paragraphs:
        if p._p.findall('.//' + qn('a:blip')):
            p.alignment = 1
            p.paragraph_format.line_spacing = 1.0
            for sh in p._p.findall('.//' + qn('wp:extent')):
                cx, cy = int(sh.get('cx')), int(sh.get('cy'))
                maxw = Inches(6.2).emu
                if cx > maxw:
                    sh.set('cx', str(maxw)); sh.set('cy', str(int(cy * maxw / cx)))
            for ext in p._p.findall('.//' + qn('a:ext')):
                pass

    # 行号（连续，全文）
    for section in d.sections:
        sectPr = section._sectPr
        ln = sectPr.find(qn('w:lnNumType'))
        if ln is None:
            ln = sectPr.makeelement(qn('w:lnNumType'), {})
            sectPr.append(ln)
        ln.set(qn('w:countBy'), '1')
        ln.set(qn('w:restart'), 'continuous')
        ln.set(qn('w:distance'), '360')

    d.save(path)
    return d


# ---------------------------------------------------------------- 阅读版（任务 #33）
NAVY, NAVY2, GREY = RGBColor(0x1F, 0x3B, 0x63), RGBColor(0x2A, 0x55, 0x8C), RGBColor(0x6B, 0x6B, 0x6B)
NAVY_HEX, ZEBRA_HEX = '1F3B63', 'EDF1F7'


def _shade(cell, hexfill):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement('w:shd')
    sh.set(qn('w:val'), 'clear'); sh.set(qn('w:color'), 'auto'); sh.set(qn('w:fill'), hexfill)
    tcPr.append(sh)


def polish_read(path):
    """阅读版：1.35 倍行距、无行号、两端对齐、标题黑体深蓝分级、表格深蓝表头 + 斑马行、
    图注 8.5 pt 灰、图居中限宽 16.3 cm。"""
    d = docx.Document(path)
    body_font = 'Times New Roman'
    head_font = 'Arial'
    cjk_body, cjk_head = ('宋体', '黑体') if ZH else (None, '黑体')

    def set_font(style, name, size, cjk=None, color=None, bold=None):
        style.font.name = name
        style.font.size = Pt(size)
        if color is not None:
            style.font.color.rgb = color
        if bold is not None:
            style.font.bold = bold
        rpr = style.element.get_or_add_rPr()
        rf = rpr.find(qn('w:rFonts'))
        if rf is None:
            rf = rpr.makeelement(qn('w:rFonts'), {}); rpr.append(rf)
        rf.set(qn('w:ascii'), name); rf.set(qn('w:hAnsi'), name)
        if cjk:
            rf.set(qn('w:eastAsia'), cjk)

    normal = d.styles['Normal']
    set_font(normal, body_font, 11.5, cjk_body)
    normal.paragraph_format.line_spacing = 1.35
    normal.paragraph_format.space_after = Pt(8)
    normal.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    for name, size, color in [('Title', 18, NAVY), ('Heading 1', 15, NAVY),
                              ('Heading 2', 13, NAVY2), ('Heading 3', 12, NAVY2)]:
        try:
            st = d.styles[name]
        except KeyError:
            continue
        set_font(st, head_font, size, cjk_head, color, True)
        st.paragraph_format.line_spacing = 1.2
        st.paragraph_format.space_before = Pt(14 if name in ('Title', 'Heading 1') else 10)
        st.paragraph_format.space_after = Pt(6)
        st.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT

    # 表格：深蓝表头（白字）+ 斑马行
    for t in d.tables:
        try:
            t.style = d.styles['Table Grid']
        except KeyError:
            pass
        for i, row in enumerate(t.rows):
            for cell in row.cells:
                if i == 0:
                    _shade(cell, NAVY_HEX)
                elif i % 2 == 0:
                    _shade(cell, ZEBRA_HEX)
                for para in cell.paragraphs:
                    para.paragraph_format.line_spacing = 1.0
                    para.paragraph_format.space_after = Pt(0)
                    para.alignment = WD_ALIGN_PARAGRAPH.LEFT
                    for r in para.runs:
                        r.font.size = Pt(9.5)
                        r.font.name = body_font
                        if i == 0:
                            r.font.bold = True
                            r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    # 图片居中、限宽 16.3 cm；图注 / 表题 8.5 pt 灰
    maxw = Cm(16.3).emu
    cap_re = re.compile(r'^\*\*(Figure|Table|图|表) ?S?\d+[.．]')
    for para in d.paragraphs:
        if para._p.findall('.//' + qn('a:blip')):
            para.alignment = WD_ALIGN_PARAGRAPH.CENTER
            para.paragraph_format.line_spacing = 1.0
            para.paragraph_format.space_after = Pt(4)
            for ext in para._p.findall('.//' + qn('wp:extent')):
                cx, cy = int(ext.get('cx')), int(ext.get('cy'))
                if cx > maxw:
                    ext.set('cx', str(maxw)); ext.set('cy', str(int(cy * maxw / cx)))
            for ext in para._p.findall('.//' + qn('a:ext')):
                cx, cy = int(ext.get('cx')), int(ext.get('cy'))
                if cx > maxw:
                    ext.set('cx', str(maxw)); ext.set('cy', str(int(cy * maxw / cx)))
            continue
        txt = para.text.strip()
        if re.match(r'^(Figure|Table|图|表) ?S?\d+[.．]', txt) or cap_re.match(txt):
            para.paragraph_format.line_spacing = 1.15
            para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            for r in para.runs:
                r.font.size = Pt(8.5)
                r.font.color.rgb = GREY

    # 阅读版不要行号
    for section in d.sections:
        sectPr = section._sectPr
        ln = sectPr.find(qn('w:lnNumType'))
        if ln is not None:
            sectPr.remove(ln)

    d.save(path)
    return d


if __name__ == '__main__':
    md = build_export_md(SRC, TMP)
    print('export md:', md)
    run_pandoc(md, OUT)
    print('pandoc ->', OUT)
    if READ:
        d = polish_read(OUT)
        print(f'  段落 {len(d.paragraphs)}，表格 {len(d.tables)}，图片 {len(d.inline_shapes)}；'
              f'阅读版：1.35 倍行距、无行号、两端对齐、深蓝分级标题、斑马表格、图注 8.5pt 灰、图限宽 16.3cm')
    else:
        d = polish(OUT)
        print(f'  段落 {len(d.paragraphs)}，表格 {len(d.tables)}，'
              f'图片 {len(d.inline_shapes)}；行号已开，双倍行距，{"宋体/" if ZH else ""}Times New Roman 12pt')
