"""Pinned finite scalar construction and ordinary-register compiler.

The source is a declared free scalar sector with supplied local services.
Its algebra and star isometries do not establish intrinsic gauge completion,
a universal solver, or autonomous control of those services.
"""


def replay(include_stars=True):
    """Run the imported compiler certificates without writing artifacts."""
    from .compiler import run as compile_cell
    from .singlet_fusion import run as compile_stars

    if type(include_stars) is not bool:
        raise ValueError("include_stars must be a boolean")
    result = {"compiler": compile_cell()}
    if include_stars:
        result["singlet_fusion"] = compile_stars()
    return result
