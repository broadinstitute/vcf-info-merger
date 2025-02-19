import logging
from pathlib import Path
from typing import Annotated, Any

import typer

from vcf_info_merger import info_merge_vcfs

app = typer.Typer()

config: dict[str, Any] = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


@app.callback(result_callback=done)
def main():
    logging.basicConfig(
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


@app.command()
def merge(
    vcf: Annotated[list[Path], typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
    chunk_size: Annotated[int, typer.Option()] = 1000000000,
) -> None:
    info_merge_vcfs(vcf_paths=vcf, out_path=out, chunk_size=chunk_size)


if __name__ == "__main__":
    app()
