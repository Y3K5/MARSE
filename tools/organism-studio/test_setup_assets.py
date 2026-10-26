"""Temporary-fixture tests; no real structural inputs modified."""

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class AssetSetupTests(unittest.TestCase):
    def test_import_and_refusal(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            target = base / "target"
            source = base / "source"
            target.mkdir()
            source.mkdir()
            shutil.copyfile(Path(__file__).with_name("setup_assets.py"), target / "setup_assets.py")
            payload = b"software fixture"
            (source / "sample.bin").write_bytes(payload)
            manifest = {
                "files": {
                    "sample.bin": {
                        "bytes": len(payload),
                        "sha256": hashlib.sha256(payload).hexdigest(),
                    }
                }
            }
            (target / "ASSET_MANIFEST.json").write_text(json.dumps(manifest))

            def run(*args):
                return subprocess.run(
                    [sys.executable, str(target / "setup_assets.py"), *args], capture_output=True
                )

            assert run().returncode != 0
            assert run("--from", str(source)).returncode == 0
            assert (target / "sample.bin").read_bytes() == payload
            (source / "sample.bin").write_bytes(b"corrupted")
            assert run("--from", str(source)).returncode != 0
            assert (target / "sample.bin").read_bytes() == payload
            assert run().returncode == 0


if __name__ == "__main__":
    unittest.main()
