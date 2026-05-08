import doctest

from hicstuff import (
    digest,
    filter,
    hicstuff,
    io,
    iteralign,
    log,
    pipeline,
    view,
)


def test_doctest():
    doctest.testmod(digest)
    doctest.testmod(filter)
    doctest.testmod(hicstuff)
    doctest.testmod(io)
    doctest.testmod(iteralign)
    doctest.testmod(log)
    doctest.testmod(pipeline)
    doctest.testmod(view)
