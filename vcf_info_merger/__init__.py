from importlib import metadata as importlib_metadata

from .merge import info_merge_vcfs


def get_version() -> str:
    try:
        return importlib_metadata.version(__name__)
    except importlib_metadata.PackageNotFoundError:
        return "unknown"


version: str = get_version()
