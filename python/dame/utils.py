import gzip


def smart_open(path, mode="r"):
    """Open a plain-text or gzip-compressed file transparently.

    Returns a text-mode file object regardless of compression, chosen by the
    ``.gz`` filename suffix. This mirrors the Rust binary, which reads gzip
    FASTQ transparently via needletail.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)
