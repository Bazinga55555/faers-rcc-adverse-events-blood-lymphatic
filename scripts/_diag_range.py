# -*- coding: utf-8 -*-
"""_diag_range.py - diagnostic only: locate the corrupt span in 2022q1.zip (no download)"""
import os
import zipfile

from paths import ROOT, DIR_RAW
zp = os.path.join(ROOT, DIR_RAW, "2022q1.zip")

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

print("corrupt entries: %d" % len(bad_entries))
for n in bad_entries:
    print("   ", n)

z = zipfile.ZipFile(zp)
spans = [(n, i.header_offset, i.compress_size) for n, i in
         ((n, z.getinfo(n)) for n in z.namelist() if not n.endswith("/"))]
z.close()

print("\nall entry offsets:")
for n, off, csz in spans:
    mark = "BAD" if n in bad_entries else "ok"
    print("   %-28s off=%12d  comp=%12d  %s" % (n, off, csz, mark))

# locate the contiguous corrupt span
first_bad_off = None
for n, off, csz in spans:
    if n in bad_entries:
        first_bad_off = off
        break

end_off = None
seen_bad = False
for n, off, csz in spans:
    if n in bad_entries:
        seen_bad = True
        continue
    if seen_bad:
        end_off = off
        break
if end_off is None:
    end_off = os.path.getsize(zp)

print("\n=== range to re-download ===")
print("   bytes %d - %d   (%.1f MB, %.1f%% of file)"
      % (first_bad_off, end_off, (end_off - first_bad_off) / 1024 / 1024,
         (end_off - first_bad_off) / os.path.getsize(zp) * 100))
