version 1.0

## Workflow: basic statistics and T2T contig detection for a genome assembly.
## Input: one assembly FASTA (plain or .gz). Output: stats table, per-contig T2T table, T2T count.

task AssemblyStatsAndT2T {
  input {
    File assembly_fasta
    Int window = 20000
    Int min_telomere_repeats = 50
    Int min_contig_length = 5000000
    Int memory_gb = 8
    Int disk_gb = 30
  }

  command <<<
    set -euo pipefail
    python3 - <<'PY'
    import gzip

    FA = "~{assembly_fasta}"
    WINDOW = ~{window}
    MIN_REP = ~{min_telomere_repeats}
    MIN_LEN = ~{min_contig_length}

    opener = gzip.open if FA.endswith(".gz") else open
    lengths = []
    t2t_rows = []

    def telo(w):
        return w.count("CCCTAA") + w.count("TTAGGG")

    def finish(name, chunks):
        seq = "".join(chunks).upper()
        lengths.append(len(seq))
        if len(seq) >= MIN_LEN:
            s = telo(seq[:WINDOW])
            e = telo(seq[-WINDOW:])
            n = seq.count("N")
            ok = s >= MIN_REP and e >= MIN_REP and n == 0
            t2t_rows.append((name, len(seq), n, s, e, "YES" if ok else "no"))

    name, chunks = None, []
    with opener(FA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    finish(name, chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        finish(name, chunks)

    lengths.sort(reverse=True)
    total = sum(lengths)
    acc, n50, l50 = 0, 0, 0
    for i, L in enumerate(lengths, start=1):
        acc += L
        if acc >= total / 2:
            n50, l50 = L, i
            break

    with open("assembly_stats.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"contigs\t{len(lengths)}\n")
        out.write(f"total_length\t{total}\n")
        out.write(f"largest_contig\t{lengths[0] if lengths else 0}\n")
        out.write(f"N50\t{n50}\n")
        out.write(f"L50\t{l50}\n")

    n_t2t = sum(1 for r in t2t_rows if r[5] == "YES")
    with open("t2t_contigs.tsv", "w") as out:
        out.write("name\tlength\tN_count\tstart_telo\tend_telo\tT2T\n")
        for r in t2t_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    with open("n_t2t.txt", "w") as out:
        out.write(str(n_t2t) + "\n")
    PY
  >>>

  output {
    File stats_tsv = "assembly_stats.tsv"
    File t2t_tsv = "t2t_contigs.tsv"
    Int n_t2t = read_int("n_t2t.txt")
  }

  runtime {
    docker: "python:3.11-slim"
    memory: "~{memory_gb} GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
  }
}

workflow AssemblyT2TReport {
  input {
    File assembly_fasta
    Int window = 20000
    Int min_telomere_repeats = 50
    Int min_contig_length = 5000000
  }

  call AssemblyStatsAndT2T {
    input:
      assembly_fasta = assembly_fasta,
      window = window,
      min_telomere_repeats = min_telomere_repeats,
      min_contig_length = min_contig_length
  }

  output {
    File stats_tsv = AssemblyStatsAndT2T.stats_tsv
    File t2t_tsv = AssemblyStatsAndT2T.t2t_tsv
    Int n_t2t = AssemblyStatsAndT2T.n_t2t
  }
}
