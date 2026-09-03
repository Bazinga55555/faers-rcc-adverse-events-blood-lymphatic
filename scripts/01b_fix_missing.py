# -*- coding: utf-8 -*-
"""
01b_fix_missing.py -- re-download the missing quarters (2017q3, 2025q3)

The earlier downloader probed remote_size unreliably while resuming with HTTP
Range requests, so some files ended up truncated (2025q3.zip was 64 MB instead
of 73 MB). This script downloads each file in a single pass and checks the size
strictly before extracting it.
Usage: python 01b_fix_missing.py
"""
import os
import io
import sys
import time
import zipfile
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from paths import ROOT, DIR_EXTRACTED, DIR_RAW
DEST = os.path.join(ROOT, DIR_RAW)
EXTR = os.path.join(ROOT, DIR_EXTRACTED)
BASE = "https://fis.fda.gov/content/Exports"

MISSING = ["2022q1"]

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
opener.addheaders = [("User-Agent", "Mozilla/5.0")]
urllib.request.install_opener(opener)


def remote_size(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    for _ in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                cr = r.headers.get("Content-Range", "")
                if "/" in cr:
                    return int(cr.split("/")[-1])
                return int(r.headers.get("Content-Length", 0) or 0)
        except Exception:
            time.sleep(3)
    return 0


def download_full(url, path, exp):
    """Download the whole file in one pass (no resumption), checking the running byte count as it goes"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for attempt in range(1, 9):
        got = 0
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(path, "wb") as f:
                while True:
                    chunk = r.read(262144)
                    if not chunk:
                        break
                    f.write(chunk)
                    got += len(chunk)
        except Exception as e:
            print("  attempt %d err=%s got=%d/%d" % (attempt, type(e).__name__, got, exp))
            time.sleep(5)
            continue
        # Verify
        actual = os.path.getsize(path)
        if actual >= exp:
            print("  complete: %d bytes" % actual)
            return True
        else:
            print("  truncated: got %d / expected %d, re-downloading" % (actual, exp))
            try:
                os.remove(path)
            except Exception:
                pass
    return False


def main():
    for q in MISSING:
        od = os.path.join(EXTR, q)
        if os.path.isdir(od) and os.listdir(od):
            print("%s already in place, skipping" % q)
            continue

        url = "%s/faers_ascii_%s.zip" % (BASE, q)
        zp = os.path.join(DEST, q + ".zip")

        print("[%s] Downloading %s ..." % (time.strftime("%H:%M:%S"), q))
        exp = remote_size(url)
        print("  remote size: %d bytes" % exp)
        if exp <= 0:
            print("  !! cannot determine the size, skipping")
            continue

        if not download_full(url, zp, exp):
            print("  !! download failed: %s" % q)
            continue

        # Extract
        try:
            if os.path.isdir(od):
                import shutil
                shutil.rmtree(od)
            os.makedirs(od, exist_ok=True)
            with zipfile.ZipFile(zp) as z:
                z.extractall(od)
            files = []
            for root, _, fs in os.walk(od):
                files.extend(fs)
            print("  extracted OK: %d files" % len(files))
            os.remove(zp)
        except Exception as e:
            print("  !! extraction failed %s: %s" % (q, repr(e)[:100]))

    print("[%s] Re-download finished" % time.strftime("%H:%M:%S"))


if __name__ == "__main__":
    main()
