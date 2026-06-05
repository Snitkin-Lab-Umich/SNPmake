#!/usr/bin/env python3
import argparse
import json

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--regions_json", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    with open(args.regions_json) as f:
        regions = json.load(f)

    with open(args.out, "w") as out:
        out.write("region\tstart\tstop\tcompleteness\tmost_common_phage\tGC\n")
        for r in regions:
            out.write(
                f"{r.get('region','')}\t{r.get('start','')}\t{r.get('stop','')}\t"
                f"{r.get('completeness','')}\t{r.get('most_common_phage','')}\t{r.get('GC','')}\n"
            )

if __name__ == "__main__":
    main()