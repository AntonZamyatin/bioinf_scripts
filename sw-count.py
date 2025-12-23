#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import pysam
import pyBigWig

# -----------------------------
# Metrics (computed per-window in one pass)
# -----------------------------

@dataclass(frozen=True)
class WindowStats:
    length: int
    a: int
    c: int
    g: int
    t: int
    n: int
    other: int

    @property
    def acgt(self) -> int:
        return self.a + self.c + self.g + self.t


def compute_window_stats(seq: str) -> WindowStats:
    a = c = g = t = n = other = 0
    for ch in seq.upper():
        if ch == "A":
            a += 1
        elif ch == "C":
            c += 1
        elif ch == "G":
            g += 1
        elif ch == "T":
            t += 1
        elif ch == "N":
            n += 1
        else:
            other += 1
    return WindowStats(length=len(seq), a=a, c=c, g=g, t=t, n=n, other=other)


MetricStatsFn = Callable[[WindowStats], float]
MetricSeqFn   = Callable[[WindowStats, str], float]

@dataclass(frozen=True)
class MetricDef:
    name: str
    fn_stats: Optional[MetricStatsFn] = None
    fn_seq: Optional[MetricSeqFn] = None

    def compute(self, st: WindowStats, seq: str) -> float:
        if self.fn_stats is not None:
            return float(self.fn_stats(st))
        if self.fn_seq is not None:
            return float(self.fn_seq(st, seq))
        raise RuntimeError(f"MetricDef {self.name} has no function")

# MetricFn = Callable[[WindowStats], float]


def metric_gc(st: WindowStats) -> float:
    denom = st.acgt
    if denom == 0:
        return float("nan")
    return (st.g + st.c) / denom


def metric_at(st: WindowStats) -> float:
    denom = st.acgt
    if denom == 0:
        return float("nan")
    return (st.a + st.t) / denom


def metric_n_frac(st: WindowStats) -> float:
    if st.length == 0:
        return float("nan")
    return st.n / st.length

import math

def tetramer_entropy_norm(st: WindowStats, seq: str, min_valid_kmers: int = 32) -> float:
    """
    Normalized Shannon entropy of 4-mer distribution in the window, in [0,1].
    Uses overlapping 4-mers, skips kmers containing non-ACGT.
    Returns NaN if too few valid 4-mers.
    """
    # 2-bit encoding for A,C,G,T
    enc = {"A": 0, "C": 1, "G": 2, "T": 3}

    counts = [0] * 256
    valid = 0

    s = seq.upper()
    n = len(s)
    if n < 4:
        return float("nan")

    code = 0
    good_run = 0  # number of consecutive valid A/C/G/T seen so far

    for ch in s:
        v = enc.get(ch)
        if v is None:
            # reset if ambiguous
            code = 0
            good_run = 0
            continue

        # rolling 2-bit code over last 4 bases
        code = ((code << 2) | v) & 0xFF
        good_run += 1
        if good_run >= 4:
            counts[code] += 1
            valid += 1

    if valid < min_valid_kmers:
        return float("nan")

    # Shannon entropy (base 2)
    inv = 1.0 / valid
    H = 0.0
    for c in counts:
        if c:
            p = c * inv
            H -= p * math.log2(p)

    # normalize by log2(256) = 8
    return H / 8.0


METRICS: Dict[str, MetricDef] = {
    "gc": MetricDef(name="gc", fn_stats=metric_gc),
    "at": MetricDef(name="at", fn_stats=metric_at),
    "n":  MetricDef(name="n",  fn_stats=metric_n_frac),
    "tetnucH": MetricDef(name="tetnucH", fn_seq=tetramer_entropy_norm)
    # add tetramer entropy below
}


# -----------------------------
# Sliding windows
# -----------------------------

@dataclass(frozen=True)
class WindowSpec:
    window: int
    step: int
    drop_last_partial: bool = True


def iter_windows(chrom_len: int, spec: WindowSpec) -> Iterable[Tuple[int, int]]:
    w, s = spec.window, spec.step
    if w <= 0 or s <= 0:
        raise ValueError("window and step must be positive integers")
    start = 0
    while start < chrom_len:
        end = start + w
        half = start + w // 2
        if end > chrom_len:
            if spec.drop_last_partial:
                break
            if half <= chrom_len:
                end = chrom_len
        yield start, end
        start += s


# -----------------------------
# File / naming
# -----------------------------

def fasta_basename_no_ext(path: str) -> str:
    name = Path(path).name
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in (".fa", ".fna", ".fasta"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


# -----------------------------
# External tools + indexing
# -----------------------------

def require_tool(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise RuntimeError(
            f"Required tool '{name}' not found in PATH. "
            f"Install it (e.g., samtools/htslib) and retry."
        )
    return exe


def run(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  {' '.join(cmd)}\n\n"
            f"STDOUT:\n{p.stdout}\n\n"
            f"STDERR:\n{p.stderr}\n"
        )


def is_bgzf_gz(path: str) -> bool:
    """
    Best-effort BGZF detection.
    BGZF is a GZIP with FEXTRA set and a 'BC' subfield in the extra area.

    This detection is sufficient for practical routing:
      - if True: treat as bgzip-compatible random-access
      - if False: treat as generic gzip (stream-only)
    """
    with open(path, "rb") as f:
        h = f.read(64)
    if len(h) < 18:
        return False
    # GZIP magic
    if h[0:2] != b"\x1f\x8b":
        return False
    flg = h[3]
    # FEXTRA flag must be set for BGZF
    if (flg & 0x04) == 0:
        return False
    # XLEN is at bytes 10-11 little-endian
    xlen = int.from_bytes(h[10:12], "little", signed=False)
    extra = h[12:12 + xlen]
    # Extra is a sequence of subfields: SI1 SI2 SLEN(2) ...data...
    i = 0
    while i + 4 <= len(extra):
        si1 = extra[i:i+1]
        si2 = extra[i+1:i+2]
        slen = int.from_bytes(extra[i+2:i+4], "little", signed=False)
        data = extra[i+4:i+4+slen]
        if si1 == b"B" and si2 == b"C" and slen == 2:
            return True
        i += 4 + slen
    return False


def ensure_fai_index(fasta_path: str) -> None:
    """
    Ensure samtools faidx index exists for fasta_path.
    Works for uncompressed FASTA and BGZF-compressed FASTA (.gz via bgzip).
    """
    require_tool("samtools")
    fai = fasta_path + ".fai"
    if os.path.exists(fai):
        return
    run(["samtools", "faidx", fasta_path])
    if not os.path.exists(fai):
        raise RuntimeError(f"samtools faidx ran but did not create index: {fai}")


def normalize_input_to_indexable_fasta(
    input_path: str,
    *,
    cache_dir: Path,
    keep_cache: bool,
) -> Tuple[str, Optional[Path]]:
    """
    Returns (indexable_fasta_path, cache_path_to_delete_or_None).

    Logic:
      - If input is not .gz: use it directly (create .fai if missing).
      - If input is .gz:
          * If BGZF: use directly (create .fai if missing).
          * Else (generic gzip): create BGZF copy in cache_dir, index that, use copy.
    """
    p = Path(input_path)
    if not p.exists():
        raise RuntimeError(f"Input file not found: {input_path}")

    if not input_path.endswith(".gz"):
        ensure_fai_index(input_path)
        return input_path, None

    # .gz: check bgzf
    if is_bgzf_gz(input_path):
        ensure_fai_index(input_path)
        return input_path, None

    # generic gzip -> make bgzip copy
    require_tool("bgzip")
    require_tool("samtools")

    cache_dir.mkdir(parents=True, exist_ok=True)
    base = fasta_basename_no_ext(input_path)
    # Put in cache_dir to avoid polluting source directory
    out_bgz = cache_dir / f"{base}.bgz.fa.gz"

    # Create BGZF: gunzip -c | bgzip -c > out_bgz
    # Use Python gzip for portability (still uses bgzip for BGZF)
    with gzip.open(input_path, "rb") as fin, open(out_bgz, "wb") as fout_raw:
        # We cannot write BGZF with Python; instead pipe into bgzip.
        # We'll stream-decompress into bgzip stdin.
        proc = subprocess.Popen(
            ["bgzip", "-c"],
            stdin=subprocess.PIPE,
            stdout=fout_raw,
            stderr=subprocess.PIPE,
        )
        assert proc.stdin is not None
        try:
            while True:
                chunk = fin.read(8 * 1024 * 1024)
                if not chunk:
                    break
                proc.stdin.write(chunk)
            proc.stdin.close()
            _, err = proc.communicate()
        except Exception:
            proc.kill()
            proc.wait()
            raise
        if proc.returncode != 0:
            raise RuntimeError(f"bgzip failed while converting gzip FASTA:\n{err.decode('utf-8', 'replace')}")

    ensure_fai_index(str(out_bgz))
    return str(out_bgz), (None if keep_cache else out_bgz)


# -----------------------------
# BigWig writing from indexed FASTA (pysam)
# -----------------------------

def write_bigwigs_indexed_fasta(
    *,
    fasta_path: str,
    outdir: Path,
    base: str,
    metric_names: List[str],
    spec: WindowSpec,
    fetch_chunk: int = 2_000_000,
    batch_windows: int = 50_000,
) -> None:
    """
    Efficient approach:
      - Get contig lengths from .fai via pysam
      - For each contig, fetch large chunks and slice windows from buffer
      - Compute WindowStats once per window, apply all metrics, write to separate BigWigs
    """
    outdir.mkdir(parents=True, exist_ok=True)

    fasta = pysam.FastaFile(fasta_path)
    try:
        chroms = [(r, fasta.get_reference_length(r)) for r in fasta.references]
        if not chroms:
            raise RuntimeError("No contigs found in FASTA.")

        # Open BigWigs
        bws: Dict[str, pyBigWig.pyBigWig] = {}
        for m in metric_names:
            out_path = outdir / f"{base}_{m}.bw"
            bw = pyBigWig.open(str(out_path), "w")
            bw.addHeader(chroms)
            bws[m] = bw

        metric_defs = {m: METRICS[m] for m in metric_names}

        try:
            for chrom, chrom_len in chroms:
                # rolling buffer for seq chunk
                buf_start = 0
                buf_end = 0
                buf_seq = ""

                # per-metric batching buffers
                starts: List[int] = []
                ends: List[int] = []
                vals: Dict[str, List[float]] = {m: [] for m in metric_names}

                def flush():
                    if not starts:
                        return
                    chroms_rep = [chrom] * len(starts)
                    for m in metric_names:
                        bws[m].addEntries(
                            chroms_rep,
                            starts,
                            ends=ends,
                            values=vals[m],
                        )
                        vals[m].clear()
                    starts.clear()
                    ends.clear()

                for start, end in iter_windows(chrom_len, spec):
                    if not (buf_start <= start and end <= buf_end):
                        buf_start = start
                        buf_end = min(chrom_len, start + fetch_chunk)
                        buf_seq = fasta.fetch(chrom, buf_start, buf_end)

                    rel_s = start - buf_start
                    rel_e = end - buf_start
                    window_seq = buf_seq[rel_s:rel_e]

                    st = compute_window_stats(window_seq)

                    starts.append(start)
                    ends.append(end)

                    for m, mdef in metric_defs.items():
                        vals[m].append(mdef.compute(st, window_seq))

                    if len(starts) >= batch_windows:
                        flush()

                flush()

        finally:
            for bw in bws.values():
                try:
                    bw.close()
                except Exception:
                    pass

    finally:
        fasta.close()


# -----------------------------
# CLI
# -----------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sw-count.py",
        description="Sliding-window genome metrics to BigWig (auto-index, supports gz via bgzip).",
    )
    p.add_argument("input_file", help="Input FASTA: .fa/.fasta/.fna or .gz")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: .)")
    p.add_argument(
        "-m", "--metrics", required=True,
        help=f"Comma-separated metric names. Available: {', '.join(sorted(METRICS))}",
    )
    p.add_argument("-w", "--window", type=int, default=1000, help="Window size (bp)")
    p.add_argument("-s", "--step", type=int, default=100, help="Step size (bp)")
    p.add_argument("--keep-last-partial", action="store_false",
                   help="Include last partial window at contig end (default: drop)")
    p.add_argument("--fetch-chunk", type=int, default=2_000_000,
                   help="Chunk size for FASTA fetch (bp) (default: 2,000,000)")
    p.add_argument("--batch-windows", type=int, default=50_000,
                   help="Batch size (windows) for BigWig addEntries (default: 50,000)")
    p.add_argument("--cache-dir", default=None,
                   help="Where to store bgzip-converted FASTA when input is generic gzip (default: temp dir)")
    p.add_argument("--keep-cache", action="store_true",
                   help="Keep bgzip-converted cached FASTA instead of deleting it")
    p.add_argument("--list-metrics", action="store_true", help="List metrics and exit")
    return p


def parse_metrics_list(raw: str) -> List[str]:
    return [x.strip() for x in raw.split(",") if x.strip()]


def main(argv: Optional[List[str]] = None) -> int:
    args = build_argparser().parse_args(argv)

    if args.list_metrics:
        for k in sorted(METRICS):
            print(k)
        return 0

    metric_names = parse_metrics_list(args.metrics)
    if not metric_names:
        print("ERROR: empty -m/--metrics list", file=sys.stderr)
        return 2
    unknown = [m for m in metric_names if m not in METRICS]
    if unknown:
        print(
            f"ERROR: Unknown metrics: {', '.join(unknown)}. "
            f"Available: {', '.join(sorted(METRICS))}",
            file=sys.stderr,
        )
        return 2

    spec = WindowSpec(
        window=int(args.window),
        step=int(args.step),
        drop_last_partial=not bool(args.keep_last_partial),
    )

    outdir = Path(args.outdir)
    base = fasta_basename_no_ext(args.input_file)

    # Choose cache dir
    if args.cache_dir:
        cache_dir = Path(args.cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        tmp_created = None
    else:
        # temp cache; can still be kept if --keep-cache
        tmp_created = tempfile.TemporaryDirectory(prefix="swcount_cache_")
        cache_dir = Path(tmp_created.name)

    # Normalize input to indexable FASTA and ensure .fai
    indexable_fasta, to_delete = normalize_input_to_indexable_fasta(
        args.input_file,
        cache_dir=cache_dir,
        keep_cache=bool(args.keep_cache),
    )

    try:
        write_bigwigs_indexed_fasta(
            fasta_path=indexable_fasta,
            outdir=outdir,
            base=base,
            metric_names=metric_names,
            spec=spec,
            fetch_chunk=int(args.fetch_chunk),
            batch_windows=int(args.batch_windows),
        )
    finally:
        # Remove cached bgzip file unless keeping
        if to_delete is not None:
            try:
                # delete .fai and .gzi if created as well
                for suf in ("", ".fai", ".gzi"):
                    p = Path(str(to_delete) + suf)
                    if p.exists():
                        p.unlink()
            except Exception:
                pass
        if tmp_created is not None:
            tmp_created.cleanup()

    for m in metric_names:
        print(outdir / f"{base}_{m}.bw")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
