"""Quiet loopback-only viewer server; no request logs or uploaded-input endpoint.

Optional private extensions: `--private DIR` (or the MARSE_PRIVATE_DIR environment variable)
mounts a folder you keep OUTSIDE this repository at /private/, read-only. The studio then reads
private/extensions.json, if present, to add links to private workspaces (extensions.js).
Without it, /private/ does not exist and nothing private is reachable or referenced. A folder
inside the repository is refused, so private data cannot be committed by accident.
"""

import argparse
import os
import posixpath
import urllib.parse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
MOUNT = "/private"
PORT = 8988


class Handler(SimpleHTTPRequestHandler):
    private_root = None
    hosts = (f"127.0.0.1:{PORT}", f"localhost:{PORT}")

    def log_message(self, *args):
        pass

    def translate_path(self, path):
        clean = urllib.parse.unquote(urllib.parse.urlsplit(path).path)
        if clean != MOUNT and not clean.startswith(MOUNT + "/"):
            return super().translate_path(path)
        missing = str(ROOT / ".no-private-mount")
        if self.private_root is None:
            return missing
        rel = posixpath.normpath("/" + clean[len(MOUNT) :]).lstrip("/")
        target = (self.private_root / rel).resolve()
        if target != self.private_root and self.private_root not in target.parents:
            return missing
        return str(target)

    def do_GET(self):
        if self.headers.get("Host", "") not in self.hosts:
            self.send_error(403)
            return
        super().do_GET()

    def do_HEAD(self):
        if self.headers.get("Host", "") not in self.hosts:
            self.send_error(403)
            return
        super().do_HEAD()

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


class ViewerServer(ThreadingHTTPServer):
    request_queue_size = 128
    daemon_threads = True


def private_root(value):
    """Resolve the optional private folder; refuse one inside this repository."""
    if not value:
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise SystemExit(f"Private folder not found: {path}")
    if path == REPO or REPO in path.parents:
        raise SystemExit(f"Keep private data outside this repository; refusing {path}")
    return path


if __name__ == "__main__":
    import json

    from setup_assets import verify

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--private",
        default=os.environ.get("MARSE_PRIVATE_DIR"),
        help="folder outside this repository to mount read-only at /private/ (MARSE_PRIVATE_DIR)",
    )
    ap.add_argument("--port", type=int, default=PORT, help="loopback port (default %(default)s)")
    args = ap.parse_args()
    errors = verify(ROOT, json.loads((ROOT / "ASSET_MANIFEST.json").read_text())["files"])
    if errors:
        raise SystemExit(
            "Runtime assets missing or changed. "
            "Run setup_assets.py --from <existing-studio-folder> first."
        )
    Handler.private_root = private_root(args.private)
    Handler.hosts = (f"127.0.0.1:{args.port}", f"localhost:{args.port}")
    mounted = " (private extensions mounted)" if Handler.private_root else ""
    print(f"MARSE viewer: http://127.0.0.1:{args.port}/unified.html{mounted}", flush=True)
    ViewerServer(("127.0.0.1", args.port), partial(Handler, directory=str(ROOT))).serve_forever()
