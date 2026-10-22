# Read generation and an independent trimming oracle for the trim suite.
import math
import os
import random
import tempfile

_dir = None


def tmp_path(name):
    global _dir
    if _dir is None:
        _dir = tempfile.mkdtemp(prefix="biotrim-test-")
    return os.path.join(_dir, name)


def _revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def make_reads(prefix, n, seed, adapter):
    """Write PREFIX_1.fq and PREFIX_2.fq: n pairs of 100-base reads from
    random fragments of 30-300 bases, reading into `adapter` (then random
    bases) past a short fragment's end; 10% of mates end in low quality, 5%
    have lowercase bases. 10% of adapters lose a base, 10% gain one and 10%
    have a substitution."""
    rng = random.Random(seed)
    tail = "".join(rng.choice("ACGT") for _ in range(200))
    with open(prefix + "_1.fq", "w") as f1, open(prefix + "_2.fq", "w") as f2:
        for i in range(n):
            frag = "".join(rng.choice("ACGT") for _ in range(rng.randint(30, 300)))
            for f, s in ((f1, frag), (f2, _revcomp(frag))):
                a = adapter
                u = rng.random()
                if u < 0.1:
                    k = rng.randrange(1, len(a))
                    a = a[:k] + a[k + 1 :]
                elif u < 0.2:
                    k = rng.randrange(1, len(a))
                    a = a[:k] + rng.choice("ACGT") + a[k:]
                elif u < 0.3:
                    k = rng.randrange(len(a))
                    a = a[:k] + rng.choice([c for c in "ACGT" if c != a[k]]) + a[k + 1 :]
                seq = (s + a + tail)[:100]
                if rng.random() < 0.05:
                    seq = seq.lower()
                q = [rng.randint(30, 40) for _ in seq]
                if rng.random() < 0.1:
                    k = rng.randint(1, 40)
                    q[-k:] = [rng.randint(2, 15) for _ in range(k)]
                f.write("@r%d insert=%d\n%s\n+\n%s\n" % (i, len(frag), seq, "".join(chr(x + 33) for x in q)))


def _read_fastq(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def _quality_cut(q, cutoff):
    if cutoff <= 0:
        return len(q)
    s, best, cut = 0, 0, len(q)
    for i in range(len(q) - 1, -1, -1):
        s += cutoff - (ord(q[i]) - 33)
        if s < 0:
            break
        if s > best:
            best, cut = s, i
    return cut


def _allowed(rate, n):
    return n * round(rate * 1000000) // 1000000


def _adapter_start(seq, n, adapter, rate, min_ov):
    """cutadapt's 3' adapter rule over seq[:n]: a full edit-distance matrix
    (free leading read bases), the leftmost full match or the longest
    partial match ending the read, whichever starts first; starts found by
    backtracking (diagonal, then deletion, then insertion)."""
    a, r = adapter.upper(), seq[:n].upper()
    m = len(a)
    if n == 0 or m == 0:
        return n
    d = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        d[i][0] = i
        for j in range(1, n + 1):
            d[i][j] = min(d[i - 1][j - 1] + (a[i - 1] != r[j - 1]), d[i - 1][j] + 1, d[i][j - 1] + 1)

    def back(i, j):
        while i > 0:
            if j > 0 and d[i][j] == d[i - 1][j - 1] + (a[i - 1] != r[j - 1]):
                i, j = i - 1, j - 1
            elif d[i][j] == d[i - 1][j] + 1:
                i -= 1
            else:
                j -= 1
        return j

    full = next((back(m, j) for j in range(1, n + 1) if d[m][j] <= _allowed(rate, m)), None)
    partial = next((back(i, n) for i in range(min(m, n), max(min_ov, 1) - 1, -1) if d[i][n] <= _allowed(rate, i)), None)
    starts = [x for x in (full, partial) if x is not None]
    return min(starts) if starts else n


def _expected(in1, in2, o):
    out, stats = [], [0, 0, 0, 0]
    for (h1, s1, q1), (h2, s2, q2) in zip(_read_fastq(in1), _read_fastq(in2)):
        mates = []
        for k, (h, s, q, a) in enumerate(((h1, s1, q1, o["adapter"]), (h2, s2, q2, o["adapter2"]))):
            qc = _quality_cut(q, o["quality"])
            cut = _adapter_start(s, qc, a, o["errorRate"], o["minOverlap"])
            stats[2 + k] += cut < qc
            mates.append((h, s[:cut], q[:cut]))
        stats[0] += 1
        if len(mates[0][1]) >= o["minLength"] and len(mates[1][1]) >= o["minLength"]:
            stats[1] += 1
            out.append((mates[0], mates[1]))
    return out, stats


def oracle_check(in1, in2, out1, out2, o):
    """'' when OUT holds exactly the pairs the oracle keeps, else why not."""
    want, _ = _expected(in1, in2, o)
    got = list(zip(_read_fastq(out1), _read_fastq(out2)))
    if len(got) != len(want):
        return "%d pairs written, expected %d" % (len(got), len(want))
    for i, (g, w) in enumerate(zip(got, want)):
        if g != w:
            return "pair %d: got %r, expected %r" % (i + 1, g, w)
    return ""


def oracle_stats(in1, in2, o):
    """(pairs, kept, first mates with adapter, second mates with adapter)"""
    return tuple(_expected(in1, in2, o)[1])


def _stats(in1, in2, o):
    st = {"basesIn": [0, 0], "qualityBases": [0, 0], "adapterBases": [0, 0],
          "basesOut": [0, 0], "adapterLengths": [{}, {}]}
    for (h1, s1, q1), (h2, s2, q2) in zip(_read_fastq(in1), _read_fastq(in2)):
        cuts = []
        for k, (s, q, a) in enumerate(((s1, q1, o["adapter"]), (s2, q2, o["adapter2"]))):
            qc = _quality_cut(q, o["quality"])
            cut = _adapter_start(s, qc, a, o["errorRate"], o["minOverlap"])
            st["basesIn"][k] += len(s)
            st["qualityBases"][k] += len(s) - qc
            st["adapterBases"][k] += qc - cut
            if cut < qc:
                h = st["adapterLengths"][k]
                h[qc - cut] = h.get(qc - cut, 0) + 1
            cuts.append(cut)
        if min(cuts) >= o["minLength"]:
            st["basesOut"][0] += cuts[0]
            st["basesOut"][1] += cuts[1]
    return st


def oracle_mate_stats(in1, in2, o, got):
    """'' when the per-mate statistics match the oracle, else the first
    difference."""
    want = _stats(in1, in2, o)
    for k in ("basesIn", "qualityBases", "adapterBases", "basesOut"):
        if list(got[k]) != want[k]:
            return "%s: got %r, want %r" % (k, list(got[k]), want[k])
    for m in (0, 1):
        g = {i: c for i, c in enumerate(got["adapterLengths"][m]) if c}
        if g != want["adapterLengths"][m]:
            return "adapterLengths[%d] differ" % m
    return ""


def check_trim_galore_json(text, mate, st):
    """'' when the JSON report has what MultiQC reads, agreeing with st."""
    import json
    try:
        d = json.loads(text)
    except ValueError as e:
        return "not JSON: %s" % e
    k = mate - 1
    checks = [
        (d["tool"], "Trim Galore"),
        (d["schema_version"], 1),
        (d["read_number"], mate),
        (d["mode"], "paired-end"),
        (d["read_processing"]["total_reads"], st["pairs"]),
        (d["read_processing"]["reads_written"], st["kept"]),
        (d["read_processing"]["reads_with_adapter"], st["adapterReads"][k]),
        (d["basepair_processing"]["total_bp_processed"], st["basesIn"][k]),
        (d["basepair_processing"]["quality_trimmed_bp"], st["qualityBases"][k]),
        (d["basepair_processing"]["total_bp_written"], st["basesOut"][k]),
        (d["pair_validation"]["pairs_removed"], st["pairs"] - st["kept"]),
        (sum(d["adapter_trimming"][0]["length_distribution"].values()), st["adapterReads"][k]),
    ]
    for i, (g, w) in enumerate(checks):
        if g != w:
            return "check %d: got %r, want %r" % (i, g, w)
    for key in ("poly_a_trimming", "poly_g_trimming", "rrbs", "parameters", "input_filenames"):
        if key not in d:
            return "missing %s" % key
    return ""


def oracle_single_stats(in1, o):
    """(reads, kept, reads with adapter, bases in, quality-trimmed bases,
    adapter bases, bases out) for single-end trimming of IN1."""
    out = [0] * 7
    for h, s, q in _read_fastq(in1):
        qc = _quality_cut(q, o["quality"])
        cut = _adapter_start(s, qc, o["adapter"], o["errorRate"], o["minOverlap"])
        out[0] += 1
        out[2] += cut < qc
        out[3] += len(s)
        out[4] += len(s) - qc
        out[5] += qc - cut
        if cut >= o["minLength"]:
            out[1] += 1
            out[6] += cut
    return tuple(out)


def check_single_json(text, st):
    import json
    d = json.loads(text)
    checks = [
        (d["mode"], "single-end"),
        (d["pair_validation"], None),
        (len(d["input_filenames"]), 1),
        (d["read_processing"]["total_reads"], st["pairs"]),
        (d["read_processing"]["reads_written"], st["kept"]),
        (d["basepair_processing"]["total_bp_written"], st["basesOut"][0]),
    ]
    for i, (g, w) in enumerate(checks):
        if g != w:
            return "check %d: got %r, want %r" % (i, g, w)
    return ""


def oracle_index(in1, in2, o):
    """The index of each pair the oracle keeps."""
    out = []
    for i, ((h1, s1, q1), (h2, s2, q2)) in enumerate(zip(_read_fastq(in1), _read_fastq(in2))):
        cuts = []
        for s, q, a in ((s1, q1, o["adapter"]), (s2, q2, o["adapter2"])):
            cuts.append(_adapter_start(s, _quality_cut(q, o["quality"]), a, o["errorRate"], o["minOverlap"]))
        if min(cuts) >= o["minLength"]:
            out.append(i)
    return out


def oracle_single_index(in1, o):
    """The index of each single-end read the oracle keeps."""
    out = []
    for i, (h, s, q) in enumerate(_read_fastq(in1)):
        if _adapter_start(s, _quality_cut(q, o["quality"]), o["adapter"], o["errorRate"], o["minOverlap"]) >= o["minLength"]:
            out.append(i)
    return out
