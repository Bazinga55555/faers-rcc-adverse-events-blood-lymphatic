# -*- coding: utf-8 -*-
"""Probe the aggregate download throughput of the FDA server across several
concurrent connections.
"""
import sys, time, threading, urllib.request

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
opener.addheaders = [('User-Agent', 'Mozilla/5.0')]
BASE = "https://fis.fda.gov/content/Exports"

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
DUR = 30
urls = ["%s/faers_ascii_202%dq%d.zip" % (BASE, 0 + (i // 4), i % 4 + 1) for i in range(N)]
res = {}


def dl(u):
    t = time.time()
    n = 0
    try:
        r = opener.open(u, timeout=60)
        while time.time() - t < DUR:
            b = r.read(65536)
            if not b:
                break
            n += len(b)
    except Exception as e:
        print("  ERR", u.split('_')[-1], e)
        n = 0
    res[u] = (n, time.time() - t)


print("Fetching %d files concurrently, %d s each ..." % (N, DUR))
ths = [threading.Thread(target=dl, args=(u,)) for u in urls]
for t in ths:
    t.start()
for t in ths:
    t.join()

tot = sum(v[0] for v in res.values())
el = max(v[1] for v in res.values())
for u, v in res.items():
    print("  %-10s %6.2f MB  %7.1f KB/s" % (u.split('_')[-1], v[0] / 1048576, v[0] / 1024 / v[1]))
print("Aggregate throughput: %.1f KB/s  (%.2f MB / %.1f s)" % (tot / 1024 / el, tot / 1048576, el))
print("Estimated time for the full 3.6 GB: %.1f hours" % (3.6 * 1024 * 1024 / (tot / 1024 / el) / 3600))
