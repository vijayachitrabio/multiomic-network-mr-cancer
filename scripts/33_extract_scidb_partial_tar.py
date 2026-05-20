#!/usr/bin/env python3

import argparse
import io
import os
import tarfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract selected files for one protein from a partial SciDB tar archive."
    )
    parser.add_argument("tar_path", help="Path to the partial .tar file")
    parser.add_argument("protein", help="Protein/project folder name inside the archive, e.g. efna1")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory (default: data/scidb/<protein>)",
    )
    parser.add_argument(
        "--want",
        nargs="*",
        default=[
            "output/cs_all2.rds",
            "output/summ_all2.assoc.txt.gz",
        ],
        help="Relative paths within the protein folder to extract",
    )
    return parser.parse_args()


def safe_extract_member(tar, member, out_path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    extracted = tar.extractfile(member)
    if extracted is None:
        return False
    with extracted, open(out_path, "wb") as handle:
        while True:
            chunk = extracted.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    return True


def main():
    args = parse_args()
    tar_path = Path(args.tar_path)
    protein = args.protein.lower()
    outdir = Path(args.outdir or f"data/scidb/{protein}")
    wanted = {f"{protein}/{rel}": rel for rel in args.want}
    found = {}

    if not tar_path.exists():
        raise SystemExit(f"Missing tar file: {tar_path}")

    # Sequential streaming mode tolerates truncated tar files as long as the
    # desired members appear before EOF.
    try:
        with tarfile.open(tar_path, mode="r|") as archive:
            for member in archive:
                if member.name not in wanted:
                    continue
                rel = wanted[member.name]
                dest = outdir / rel
                ok = safe_extract_member(archive, member, dest)
                if ok:
                    found[member.name] = dest
                if len(found) == len(wanted):
                    break
    except (tarfile.TarError, OSError, EOFError) as exc:
        print(f"Stopped while scanning partial tar: {exc}")

    print("Found files:")
    for name, dest in found.items():
        print(f"- {name} -> {dest}")

    missing = [name for name in wanted if name not in found]
    if missing:
        print("Still missing:")
        for name in missing:
            print(f"- {name}")
        raise SystemExit(2)


if __name__ == "__main__":
    main()
