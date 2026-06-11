import gzip

from dame.utils import smart_open


def test_smart_open_plain(tmp_path):
    p = tmp_path / "plain.txt"
    p.write_text("hello\nworld\n")
    with smart_open(str(p)) as f:
        assert f.read() == "hello\nworld\n"


def test_smart_open_gzip(tmp_path):
    p = tmp_path / "data.txt.gz"
    with gzip.open(str(p), "wt") as f:
        f.write("hello\nworld\n")
    # Returns text (not bytes) for a .gz path
    with smart_open(str(p)) as f:
        assert f.read() == "hello\nworld\n"
