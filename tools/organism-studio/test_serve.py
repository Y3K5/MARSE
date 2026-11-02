"""Temporary-fixture tests for the optional private mount in serve.py; nothing outside a temporary
folder is read."""

import os
import tempfile
import unittest
from pathlib import Path

import serve


def handler(private=None):
    h = serve.Handler.__new__(serve.Handler)
    h.directory = str(serve.ROOT)
    h.private_root = private
    return h


def refused(value):
    try:
        serve.private_root(value)
    except SystemExit:
        return True
    return False


class PrivateMountTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.private = Path(self.tmp.name).resolve() / "private-data"
        (self.private / "pieces").mkdir(parents=True)
        (self.private / "extensions.json").write_text("{}")
        (self.private / "pieces" / "studio.html").write_text("<p>private</p>")

    def tearDown(self):
        self.tmp.cleanup()

    def test_no_mount_means_no_private_paths(self):
        for path in ("/private", "/private/extensions.json", "/private/../serve.py"):
            assert not Path(handler().translate_path(path)).exists(), path

    def test_mount_serves_files_inside_the_private_folder(self):
        h = handler(self.private)
        assert (
            Path(h.translate_path("/private/extensions.json")) == self.private / "extensions.json"
        )
        piece = h.translate_path("/private/pieces/studio.html?source=2OMF")
        assert Path(piece) == self.private / "pieces" / "studio.html"

    def test_traversal_and_escaping_symlinks_stay_inside(self):
        outside = Path(self.tmp.name).resolve() / "outside.txt"
        outside.write_text("not for the studio")
        os.symlink(outside, self.private / "escape.txt")
        h = handler(self.private)
        for path in (
            "/private/../outside.txt",
            "/private/%2e%2e/outside.txt",
            "/private/escape.txt",
        ):
            assert Path(h.translate_path(path)).resolve() != outside, path

    def test_public_paths_are_unchanged(self):
        public = handler(self.private).translate_path("/unified.html")
        assert Path(public) == serve.ROOT / "unified.html"

    def test_private_folder_inside_the_repository_is_refused(self):
        assert refused(str(serve.ROOT))
        assert refused(str(self.private / "missing"))
        assert serve.private_root(str(self.private)) == self.private
        assert serve.private_root(None) is None


if __name__ == "__main__":
    unittest.main()
