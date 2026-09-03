# -*- coding: utf-8 -*-
"""
FAERS / AEMS quarterly data downloader (8 parallel workers, resumable,
--------------------------------------------------------------
Data window : 2004Q1 - 2026Q2  (90 quarters)
URL pattern :
    2004Q1-2012Q3 -> https://fis.fda.gov/content/Exports/aers_ascii_YYYYqN.zip
    2012Q4-2026Q2 -> https://fis.fda.gov/content/Exports/faers_ascii_YYYYqN.zip
Notes:
    * The machine must connect directly (no proxy); routing through a proxy returns HTTP=000
    * The server throttles one connection to ~120-460 KB/s, but separate files can be
    * fetched in parallel, so 8 workers reach ~2 MB/s in aggregate
"""
import os, io, sys, time, zipfile, threading
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

from paths import ROOT, DIR_EXTRACTED, DIR_RAW, DIR_SCRIPTS
DEST   = os.path.join(ROOT, DIR_RAW)
EXTR   = os.path.join(ROOT, DIR_EXTRACTED)
LOG    = os.path.join(ROOT, DIR_SCRIPTS, "download.log")
BASE   = "https://fis.fda.gov/content/Exports"
END    = (2026, 2)
WORKERS = 8

os.makedirs(DEST, exist_ok=True)
os.makedirs(EXTR, exist_ok=True)

# Disable any system proxy, otherwise the connection fails
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
opener.addheaders = [("User-Agent", "Mozilla/5.0")]
urllib.request.install_opener(opener)

_lock = threading.Lock()


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with _lock:
        with io.open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def quarters():
    out = []
    for y in range(2004, 2027):
        for q in range(1, 5):
            if (y, q) > END:
                break
            pre = "aers" if (y < 2012 or (y == 2012 and q <= 3)) else "faers"
            out.append((pre, "%dq%d" % (y, q)))
    return out


def remote_size(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            cr = r.headers.get("Content-Range", "")
            if "/" in cr:
                return int(cr.split("/")[-1])
            return int(r.headers.get("Content-Length", 0) or 0)
    except Exception:
        return 0


def fetch(prefix, q, tries=8):
    url = "%s/%s_ascii_%s.zip" % (BASE, prefix, q)
    zp  = os.path.join(DEST, q + ".zip")
    od  = os.path.join(EXTR, q)

    if os.path.isdir(od) and os.listdir(od):
        return "SKIP", 0

    exp = remote_size(url)
    if exp <= 0:
        return "NOSIZE", 0

    for t in range(1, tries + 1):
        have = os.path.getsize(zp) if os.path.exists(zp) else 0
        if have >= exp:
            break
        try:
            req = urllib.request.Request(url)
            if have > 0:
                req.add_header("Range", "bytes=%d-" % have)
            with urllib.request.urlopen(req, timeout=300) as r, open(zp, "ab") as f:
                while True:
                    chunk = r.read(262144)
                    if not chunk:
                        break
                    f.write(chunk)
        except Exception as e:
            log("  retry %s t=%d err=%s" % (q, t, type(e).__name__))
            time.sleep(4)

    have = os.path.getsize(zp) if os.path.exists(zp) else 0
    if have < exp:
        return "INCOMPLETE", have

    # Extract
    try:
        os.makedirs(od, exist_ok=True)
        with zipfile.ZipFile(zp) as z:
            z.extractall(od)
        n = len(os.listdir(od))
    except Exception as e:
        return "UNZIPFAIL", have
    if n:
        try:
            os.remove(zp)
        except Exception:
            pass
    return ("OK", have) if n else ("EMPTY", have)


if __name__ == "__main__":
    log("=" * 52)
    log("Downloading 2004Q1-2026Q2 with %d parallel workers" % WORKERS)
    stat = {}
    t0 = time.time()
    todo = [(pre, q) for pre, q in quarters()
            if not (os.path.isdir(os.path.join(EXTR, q)) and os.listdir(os.path.join(EXTR, q)))]
    log("%d quarters to download out of %d" % (len(todo), len(quarters())))

    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(fetch, pre, q): q for pre, q in todo}
        for fu in as_completed(futs):
            q = futs[fu]
            try:
                s, sz = fu.result()
            except Exception as e:
                s, sz = "ERROR", 0
                log("EXC %s %s" % (q, e))
            stat[s] = stat.get(s, 0) + 1
            done += 1
            if s in ("OK", "SKIP"):
                log("%-6s %s   %6.1f MB   [%d/%d]   %.0f min elapsed"
                    % (s, q, sz / 1048576.0, done, len(todo), (time.time() - t0) / 60))
            else:
                log("%-6s %s   %d bytes" % (s, q, sz))

    log("Done: %s   total %.0f min" % (stat, (time.time() - t0) / 60))
    log("=" * 52)
