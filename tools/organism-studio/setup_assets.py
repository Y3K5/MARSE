"""Verify or import exact local studio assets. No network access or installation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

ROOT = Path(__file__).resolve().parent

def verify(root, files):
    errors = []
    for name, expected in files.items():
        rel = Path(name)
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError("Unsafe asset path")
        path = root / rel
        if path.is_symlink() or not path.is_file():
            errors.append(name + ": missing or symlink")
        elif path.stat().st_size != expected["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != expected["sha256"]:
            errors.append(name + ": checksum/size mismatch")
    return errors

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from", dest="source", type=Path, help="Existing studio folder containing vendor, structures and integration")
    args = parser.parse_args()
    files = json.loads((ROOT / "ASSET_MANIFEST.json").read_text())["files"]
    source = args.source.resolve() if args.source else ROOT
    errors = verify(source, files)
    if errors:
        raise SystemExit("Asset verification failed; nothing copied:\n" + "\n".join(errors))
    if args.source and source != ROOT:
        # Validate every source before any destination changes, then atomically replace each file.
        for name in files:
            dest = ROOT / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as tmp:
                temp = Path(tmp.name)
            try:
                shutil.copyfile(source / name, temp)
                if hashlib.sha256(temp.read_bytes()).hexdigest() != files[name]["sha256"]:
                    raise ValueError("Source changed during copy: " + name)
                temp.replace(dest)
            finally:
                temp.unlink(missing_ok=True)
    errors = verify(ROOT, files)
    if errors:
        raise SystemExit("Destination verification failed: " + "; ".join(errors))
    print(str(len(files)) + " local assets verified. No network requests made.")

if __name__ == "__main__":
    main()
