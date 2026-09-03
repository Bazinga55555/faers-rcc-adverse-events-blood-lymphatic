#!/bin/bash
# ============================================================
# FAERS / AEMS 季度数据下载（断点续传版）
#   2004Q1 - 2012Q3 : aers_ascii_YYYYqN.zip   (35 个)
#   2012Q4 - 2026Q2 : faers_ascii_YYYYqN.zip  (55 个)
#   - 本机必须直连：--noproxy '*'（走代理会 HTTP=000）
#   - 服务端限速约 190 KB/s，并发无效，故串行 + curl -C - 断点续传
# ============================================================
set -u

BASE="https://fis.fda.gov/content/Exports"
ROOT="/f/2026-2027-1科研/四类机制 × 血液和淋巴系统疾病-workbuddy"
DEST="$ROOT/01_原始数据"
EXTRACT="$ROOT/02_解压数据"
LOG="$ROOT/04_分析脚本/download.log"

mkdir -p "$DEST" "$EXTRACT"

quarters () {
    for ((y=2004; y<=2026; y++)); do
        for ((q=1; q<=4; q++)); do
            [ "$y" -eq 2026 ] && [ "$q" -gt 2 ] && break
            if [ "$y" -lt 2012 ] || { [ "$y" -eq 2012 ] && [ "$q" -le 3 ]; }; then
                echo "aers ${y}q${q}"
            else
                echo "faers ${y}q${q}"
            fi
        done
    done
}

log(){ echo "[$(date '+%H:%M:%S')] $*" >> "$LOG"; }

log "========== 开始下载 2004Q1-2026Q2 =========="

n_ok=0; n_fail=0
while read -r prefix q; do
    [ -z "${prefix:-}" ] && continue
    url="${BASE}/${prefix}_ascii_${q}.zip"
    zipfile="$DEST/${q}.zip"
    outdir="$EXTRACT/${q}"

    # 已解压成功 -> 跳过
    if [ -d "$outdir" ] && [ -n "$(ls -A "$outdir" 2>/dev/null)" ]; then
        log "SKIP  $q (已解压)"
        n_ok=$((n_ok+1)); continue
    fi

    # 断点续传重试最多 8 次
    ok=0
    for attempt in 1 2 3 4 5 6 7 8; do
        code=$(curl -s --noproxy '*' --connect-timeout 30 --max-time 1800 \
                    -C - -o "$zipfile" -w "%{http_code}" "$url" 2>/dev/null)
        sz=$(stat -c %s "$zipfile" 2>/dev/null || echo 0)

        # 拿到期望大小（只在首轮探测一次）
        if [ -z "${exp:-}" ] || [ "$attempt" -gt 1 ]; then
            exp=$(curl -s --noproxy '*' --max-time 40 -r 0-999999999 \
                  -D - -o /dev/null "$url" 2>/dev/null \
                  | grep -i '^content-range' | tr -d '\r' | sed 's#.*/##')
        fi

        if [ "$code" = "200" ] || [ "$code" = "206" ]; then
            if [ -n "${exp:-}" ] && [ "$exp" -gt 0 ] && [ "$sz" -ge "$exp" ]; then ok=1; break; fi
            if [ -z "${exp:-}" ] || [ "$exp" = "0" ]; then
                # 拿不到大小就按 zip 尾签判断
                if tail -c 22 "$zipfile" 2>/dev/null | grep -q "PK"; then ok=1; break; fi
            fi
        fi
        log "RETRY $q  try=$attempt  http=$code  got=$sz  want=${exp:-?}"
        sleep 5
    done

    if [ "$ok" -ne 1 ]; then
        log "FAIL  $q"
        n_fail=$((n_fail+1)); unset exp; continue
    fi

    mkdir -p "$outdir"
    if unzip -oq "$zipfile" -d "$outdir" 2>/dev/null && [ -n "$(ls -A "$outdir" 2>/dev/null)" ]; then
        log "OK    $q  $((sz/1048576))MB  files=$(ls -1 "$outdir" | wc -l)"
        rm -f "$zipfile"
        n_ok=$((n_ok+1))
    else
        log "FAIL  $q 解压失败"
        rm -f "$zipfile"; rm -rf "$outdir"
        n_fail=$((n_fail+1))
    fi
    unset exp
done < <(quarters)

log "========== 完成: OK=$n_ok  FAIL=$n_fail =========="
echo "DONE OK=$n_ok FAIL=$n_fail" >> "$LOG"
