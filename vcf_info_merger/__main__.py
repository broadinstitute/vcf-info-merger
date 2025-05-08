import logging
import sys
from pathlib import Path
from typing import Annotated, Any

import typer

from vcf_info_merger import info_merge_vcfs

app = typer.Typer()

config: dict[str, Any] = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


def set_up_gcp_friendly_logging(level: int = logging.INFO) -> None:
    """
    Configure logging so that logs are routed to stdout/stderr based on severity,
    for compatibility with Google Cloud Logging, and are prepended by timestamps.

    :param level: log level to set
    """

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()

    # formatter for all log levels
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
    )

    # handler for DEBUG and INFO → stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.addFilter(lambda x: x.levelno < logging.WARNING)
    stdout_handler.setFormatter(formatter)

    # handler for WARNING and above → stderr
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)

    logger.addHandler(stdout_handler)
    logger.addHandler(stderr_handler)


@app.callback(result_callback=done)
def main():
    set_up_gcp_friendly_logging()


@app.command()
def merge(
    vcf: Annotated[list[Path], typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
    chunk_size: Annotated[int, typer.Option()] = 1000000000,
) -> None:
    info_merge_vcfs(vcf_paths=vcf, out_path=out, chunk_size=chunk_size)


if __name__ == "__main__":
    app()
