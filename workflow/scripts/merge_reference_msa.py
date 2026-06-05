import argparse
import sys
from Bio import SeqIO

DEFAULT_VALID = set("ACGTN-")

def read_ref_seq(ref_fasta, ref_id):
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        if rec.id == ref_id:
            return str(rec.seq).upper()
    raise SystemExit(f"[merge] ERROR: ref_id '{ref_id}' not found in ref_fasta: {ref_fasta}")

def load_fasta(path):
    """Return (seq_dict, order_list). Uses record.id (first token in header)."""
    d = {}
    order = []
    for rec in SeqIO.parse(path, "fasta"):
        rid = rec.id
        if rid in d:
            raise SystemExit(f"[merge] ERROR: duplicate record ID '{rid}' in {path}")
        d[rid] = str(rec.seq).upper()
        order.append(rid)
    return d, order

def validate_seq(seq, want_len, where, valid_chars):
    if len(seq) != want_len:
        raise SystemExit(f"[merge] ERROR: length mismatch for {where}: got {len(seq)} expected {want_len}")
    bad = sorted(set(seq) - valid_chars)
    if bad:
        raise SystemExit(f"[merge] ERROR: invalid characters for {where}: {bad}")

def write_fasta(out_path, records, wrap=0):
    with open(out_path, "w") as out:
        for rid, seq in records:
            out.write(f">{rid}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i+wrap] + "\n")
            else:
                out.write(seq + "\n")

def main():
    ap = argparse.ArgumentParser(
        description="Merge new reference-coordinate sequences into an existing reference-coordinate MSA with strict checks."
    )
    ap.add_argument("--existing_msa", required=True)
    ap.add_argument("--ref_fasta", required=True)
    ap.add_argument("--ref_id", required=True, help="Reference record ID (must be present in both ref_fasta and existing_msa).")
    ap.add_argument("--add", nargs="+", required=True, help="One or more FASTA files containing sequences to add.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--duplicate_policy", choices=["error", "skip", "overwrite"], default="error")
    ap.add_argument("--wrap", type=int, default=0)
    ap.add_argument("--valid_chars", default="ACGTN-",
                    help="Allowed characters (default ACGTN-). Example to allow ambiguity: ACGTN-WSMKRYBDHV")
    args = ap.parse_args()

    valid_chars = set(args.valid_chars.upper())

    ref_seq = read_ref_seq(args.ref_fasta, args.ref_id)
    ref_len = len(ref_seq)

    # Load and validate existing MSA
    msa, msa_order = load_fasta(args.existing_msa)

    if args.ref_id not in msa:
        raise SystemExit(f"[merge] ERROR: ref_id '{args.ref_id}' not found in existing MSA: {args.existing_msa}")

    if msa[args.ref_id] != ref_seq:
        raise SystemExit(
            "[merge] ERROR: reference sequence mismatch between existing MSA and ref_fasta.\n"
            "This indicates different reference/coordinates. Refusing to merge."
        )

    for rid in msa_order:
        validate_seq(msa[rid], ref_len, f"{args.existing_msa}:{rid}", valid_chars)

    # Ensure reference is first in output order
    msa_order_no_ref = [rid for rid in msa_order if rid != args.ref_id]
    out_order = [args.ref_id] + msa_order_no_ref

    added = 0
    skipped = 0
    overwritten = 0

    # Add new sequences
    for fa in args.add:
        add_dict, add_order = load_fasta(fa)
        for rid in add_order:
            seq = add_dict[rid]

            # allow reference record in add files, but verify and skip
            if rid == args.ref_id:
                if seq != ref_seq:
                    raise SystemExit(f"[merge] ERROR: ref record in {fa} does not match ref_fasta ({args.ref_id})")
                continue

            validate_seq(seq, ref_len, f"{fa}:{rid}", valid_chars)

            if rid in msa:
                if args.duplicate_policy == "error":
                    raise SystemExit(f"[merge] ERROR: sample '{rid}' already exists in existing MSA (duplicate_policy=error)")
                elif args.duplicate_policy == "skip":
                    skipped += 1
                    continue
                elif args.duplicate_policy == "overwrite":
                    msa[rid] = seq
                    overwritten += 1
                    continue

            msa[rid] = seq
            out_order.append(rid)
            added += 1

    # Write output
    records = [(rid, msa[rid]) for rid in out_order]
    write_fasta(args.out, records, wrap=args.wrap)

    sys.stderr.write(
        f"[merge] existing={len(msa_order)} added={added} overwritten={overwritten} skipped={skipped} "
        f"out={len(records)} ref_len={ref_len}\n"
    )

if __name__ == "__main__":
    main()