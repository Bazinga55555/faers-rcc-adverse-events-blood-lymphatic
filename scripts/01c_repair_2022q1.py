# -*- coding: utf-8 -*-
"""
01c_repair_2022q1.py -- precisely repair the contiguous corrupt region of 2022q1.zip

Background:
  After 2022q1.zip was downloaded in full (67,926,472 bytes), unzipping it
  raised Bad CRC-32. The damage was traced to one contiguous 49 MB block:
      bytes 493,513 - 51,930,575 (about 49 MB)
  which covers the data area of ASCII/DEMO22Q1.txt and the local file headers
  of the entries that follow it (DRUG22Q1.pdf / DRUG22Q1.txt / INDI22Q1.pdf /
  INDI22Q1.txt / OUTC22Q1.pdf). Integrity resumes only at ASCII/OUTC22Q1.txt (data_start 51,930,624).

  The parts already confirmed intact (17 MB, including the core tables REAC/THER/RPSR/OUTC) are left alone.

Strategy:
  1. Re-download the damaged range with HTTP Range requests, 4 MB at a time
  2. Retry each chunk on its own (a failure then costs 4 MB, not 49 MB)
  3. Once everything is written back, re-check the CRC of every entry with zipfile
  4. Extract only if verification passes

Usage:
  python 01c_repair_2022q1.py
"""
import io
import os
import sys
import time
import zipfile
import struct
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from paths import ROOT, DIR_EXTRACTED, DIR_RAW
DEST  = os.path.join(ROOT, DIR_RAW)
EXTR  = os.path.join(ROOT, DIR_EXTRACTED)
BASE  = "https://fis.fda.gov/content/Exports"
QUARTER = "2022q1"
CHUNK = 4 * 1024 * 1024      # 4 MB per chunk
MAX_RETRY_PER_CHUNK = 6

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


def bad_ranges(zp):
    """Return (list of byte ranges to repair, names of all bad entries)"""
    z = zipfile.ZipFile(zp)
    bad_entries = []
    for n in z.namelist():
        if n.endswith("/"):
            continue
        try:
            with z.open(n) as src:
                while src.read(1 << 20):
                    pass
        except Exception:
            bad_entries.append(n)
    z.close()

    if not bad_entries:
        return [], [], []

    # Work out the data range of every bad entry, plus entries that cannot be located because their local header is damaged
    z = zipfile.ZipFile(zp)
    spans = []
    for n in z.namelist():
        if n.endswith("/"):
            continue
        i = z.getinfo(n)
        spans.append((n, i.header_offset, i.CRC, i.compress_size))
    z.close()

    # Start of the damaged range = header_offset of the first bad entry
    # End of the damaged range  = header_offset of the first later entry that still opens cleanly
    first_bad_off = None
    for n, off, crc, csz in spans:
        if n in bad_entries:
            first_bad_off = off
            break
    if first_bad_off is None:
        # The local header of the bad entry is damaged too, so infer it from where the previous entry ends
        for idx, (n, off, crc, csz) in enumerate(spans):
            if n in bad_entries:
                first_bad_off = spans[idx - 1][1] if idx > 0 else 0
                break

    end_off = None
    seen_bad = False
    for n, off, crc, csz in spans:
        if n in bad_entries:
            seen_bad = True
            continue
        if seen_bad:
            end_off = off
            break
    if end_off is None:
        end_off = os.path.getsize(zp)

    return [(first_bad_off, end_off)], bad_entries, spans


def download_range(url, start, end, path, label=""):
    """Download the byte range [start, end) and write it into path at the matching offset; retries on failure"""
    want = end - start
    for attempt in range(1, MAX_RETRY_PER_CHUNK + 1):
        got = 0
        try:
            req = urllib.request.Request(
                url, headers={"Range": "bytes=%d-%d" % (start, end - 1)})
            with urllib.request.urlopen(req, timeout=180) as r, open(path, "r+b") as f:
                f.seek(start)
                while True:
                    chunk = r.read(262144)
                    if not chunk:
                        break
                    f.write(chunk)
                    got += len(chunk)
            if got == want:
                return True, got
            print("    %s attempt %d short read got=%d/%d" % (label, attempt, got, want))
        except Exception as e:
            print("    %s attempt %d err=%s got=%d/%d"
                  % (label, attempt, type(e).__name__, got, want))
        time.sleep(3)
    return False, got


def verify(zp):
    """Return the list of bad entries"""
    z = zipfile.ZipFile(zp)
    bad = []
    for n in z.namelist():
        if n.endswith("/"):
            continue
        try:
            with z.open(n) as src:
                while src.read(1 << 20):
                    pass
        except Exception as e:
            bad.append((n, type(e).__name__))
    z.close()
    return bad


def main():
    zp = os.path.join(DEST, QUARTER + ".zip")
    url = "%s/faers_ascii_%s.zip" % (BASE, QUARTER)

    if not os.path.exists(zp):
        print("!! %s not found" % zp)
        return

    exp = remote_size(url)
    actual = os.path.getsize(zp)
    print("[%s] local %d bytes, remote %d bytes" % (time.strftime("%H:%M:%S"), actual, exp))
    if actual != exp:
        print("!! Size mismatch: the whole file has to be re-downloaded (this script only repairs part of it)")
        return

    bad = verify(zp)
    if not bad:
        print("zip is intact, nothing to repair")
    else:
        print("bad entries: %d" % len(bad))
        for n, e in bad:
            print("   %-32s %s" % (n, e))

        ranges, bad_entries, spans = bad_ranges(zp)
        print("\nRanges to re-download:")
        total = 0
        for s, e in ranges:
            print("   bytes %d - %d  (%.1f MB)" % (s, e, (e - s) / 1024 / 1024))
            total += e - s
        print("   total %.1f MB (%.1f%% of the file)" % (total / 1024 / 1024, total / actual * 100))

        # Re-download chunk by chunk
        print("\nRe-downloading in chunks of %d MB ..." % (CHUNK // 1024 // 1024))
        t0 = time.time()
        for s, e in ranges:
            pos = s
            idx = 0
            while pos < e:
                ce = min(pos + CHUNK, e)
                label = "[%d/%d]" % (idx + 1, -(-(e - s) // CHUNK))
                ok, got = download_range(url, pos, ce, zp, label)
                if not ok:
                    print("  !! chunk %s still failing after %d retries, giving up" % (label, MAX_RETRY_PER_CHUNK))
                    return
                done = (pos - s + (ce - pos)) / (e - s) * 100
                print("    %s %d-%d  ok (%.1f%%, %.1f min elapsed)"
                      % (label, pos, ce, done, (time.time() - t0) / 60), flush=True)
                pos = ce
                idx += 1

        print("\nRe-download finished in %.1f min, verifying ..." % ((time.time() - t0) / 60))
        bad2 = verify(zp)
        if bad2:
            print("!! These entries are still corrupt:")
            for n, e in bad2:
                print("   %-32s %s" % (n, e))
            return
        print("All CRC checks passed")

    # Extract
    od = os.path.join(EXTR, QUARTER)
    print("\nExtracting to %s ..." % od)
    if os.path.isdir(od):
        import shutil
        shutil.rmtree(od)
    os.makedirs(od, exist_ok=True)
    with zipfile.ZipFile(zp) as z:
        z.extractall(od)
    files = []
    for root, _, fs in os.walk(od):
        files.extend(fs)
    print("Extracted OK: %d files" % len(files))
    for f in sorted(files):
        if f.endswith(".txt"):
            fp = os.path.join(od, "ASCII", f)
            if os.path.exists(fp):
                print("   %-24s %12d bytes" % (f, os.path.getsize(fp)))
    print("\n[%s] Done" % time.strftime("%H:%M:%S"))


if __name__ == "__main__":
    main()
