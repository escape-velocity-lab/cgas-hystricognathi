"""17a 辅助：运行自编的 Windows 版 HyPhy 2.5.103（32 位 MinGW，见汇报「新坑」）。

这个编译版算完、写完 JSON 后进程会挂在退出阶段（OpenMP 线程清理），所以：
  监视输出 JSON —— 文件出现且 20 s 内大小不变、能被 json 解析 → 判定完成 → 结束进程。
用法: python 14_t17_hyphy_run.py <json 输出> <日志> <hyphy 参数…>
  例：python 14_t17_hyphy_run.py out.json out.log fel --alignment a.fasta --tree t.nwk --output out.json CPU=4
"""
import sys, os, time, json, subprocess

HY = os.environ.get('HYPHY', 'hyphy')
OUT, LOG = sys.argv[1], sys.argv[2]
args = sys.argv[3:]
env = dict(os.environ, HYPHY_PATH=os.environ.get('HYPHY_PATH', ''))
if os.path.exists(OUT):
    os.remove(OUT)
t0 = time.time()
p = subprocess.Popen([HY] + args, stdout=open(LOG, 'w', encoding='utf-8', errors='replace'), stderr=subprocess.STDOUT,
                     stdin=subprocess.DEVNULL, env=env)
last, stable = -1, 0
while True:
    time.sleep(5)
    if p.poll() is not None and not os.path.exists(OUT):
        print(f'HyPhy 退出（rc={p.returncode}）但没有输出 JSON，见 {LOG}', flush=True)
        sys.exit(1)
    if os.path.exists(OUT):
        sz = os.path.getsize(OUT)
        stable = stable + 5 if sz == last else 0
        last = sz
        if stable >= 20:
            try:
                json.load(open(OUT))
                break
            except Exception:
                stable = 0
if p.poll() is None:
    p.kill()
print(f'完成：{OUT}（{time.time() - t0:.0f} s）', flush=True)
