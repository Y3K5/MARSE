"""Quiet loopback-only viewer server; no request logs or uploaded-input endpoint."""
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from pathlib import Path
from functools import partial
ROOT=Path(__file__).resolve().parent
class Handler(SimpleHTTPRequestHandler):
    def log_message(self,*args):pass
    def do_GET(self):
        if self.headers.get('Host','')not in ['127.0.0.1:8988','localhost:8988']:
            self.send_error(403);return
        super().do_GET()
    def do_HEAD(self):
        if self.headers.get('Host','')not in ['127.0.0.1:8988','localhost:8988']:
            self.send_error(403);return
        super().do_HEAD()
    def end_headers(self):
        self.send_header('Cache-Control','no-store')
        super().end_headers()
class ViewerServer(ThreadingHTTPServer):
    request_queue_size=128
    daemon_threads=True
if __name__=='__main__':
    import json
    from setup_assets import verify
    errors=verify(ROOT,json.loads((ROOT/'ASSET_MANIFEST.json').read_text())['files'])
    if errors: raise SystemExit('Runtime assets missing or changed. Run setup_assets.py --from <existing-studio-folder> first.')
    print('MARSE viewer: http://127.0.0.1:8988/unified.html',flush=True)
    ViewerServer(('127.0.0.1',8988),partial(Handler,directory=str(ROOT))).serve_forever()
