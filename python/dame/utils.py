import gzip


def smart_open(path, mode='r'):
    """Open a plain-text or gzip-compressed file transparently.

    Always returns a text-mode file object regardless of compression.
    """
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)
