import argparse

def read_positions(path):
    pos = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            pos.add(int(line))
    return pos

def main():
    p = argparse.ArgumentParser(description="Replace bases at unmapped positions with '-' in a consensus FASTA.")
    p.add_argument("--fasta", required=True, help="Input consensus FASTA (single sequence, full-length)")
    p.add_argument("--unmapped", required=True, help="Unmapped positions file (1-based positions, one per line)")
    p.add_argument("--out", required=True, help="Output FASTA")
    args = p.parse_args()

    unmapped = read_positions(args.unmapped)

    header = None
    seq_chunks = []
    with open(args.fasta) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
            elif line:
                seq_chunks.append(line.strip())

    if header is None:
        raise ValueError(f"No FASTA header found in {args.fasta}")

    seq = list("".join(seq_chunks))
    n = len(seq)

    # Replace positions with '-' (1-based -> 0-based index)
    for p1 in unmapped:
        i = p1 - 1
        if 0 <= i < n:
            seq[i] = "-"
        # if position outside length, ignore

    seq = "".join(seq)

    with open(args.out, "w") as out:
        out.write(header + "\n")
        # wrap 60 chars/line
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

if __name__ == "__main__":
    main()