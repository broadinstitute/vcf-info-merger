import concurrent.futures
import logging
import os
import re
import subprocess
import tempfile
from csv import QUOTE_NONE
from glob import glob
from pathlib import Path
from typing import Any, Hashable, Type

import bgzip
import pandas as pd

from vcf_info_merger.types import Chunk, PanderaBaseSchema, TypedDataFrame, Vcf

logger = logging.getLogger(__name__)


def info_merge_vcfs(
    vcf_paths: list[Path | str], out_path: Path | str, chunk_size: int = 1000000000
) -> None:
    """
    Merge INFO fields from multiple VCF files (for the single sample) into a single VCF
    file.

    :param vcf_paths: a list of paths to VCF files
    :param out_path: the path to write the merged VCF file
    :param chunk_size: number of lines in the VCFs to process at a time
    """

    logger.info(f"Merging {len(vcf_paths)} VCFs into {out_path}")

    header_lines = get_header_lines(vcf_paths)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # default chunk size of 1000000000 will just split on chroms, but can be reduced
        # if memory usage is a concern
        chunks = make_chunks(header_lines, chunk_size)
        chunk_vcfs(vcf_paths, tmp_dir, chunks)
        process_chunks(header_lines, out_path, tmp_dir, chunks, len(vcf_paths))

    logger.info("Done merging")


def get_header_lines(vcf_paths: list[Path | str]) -> list[str]:
    """
    Read all of the header lines in a list of VCF files and return their union,
    retaining their order.

    :param vcf_paths: a list of paths to VCF files
    :return: a list of the distinct header lines in their original order
    """

    header_lines = []
    col_header_line = None

    for path in vcf_paths:
        # don't read more lines than necessary to get the entire header
        break_next_time = False

        with open(path, "rb") as raw:
            logger.info(f"Reading {os.path.basename(path)} header")
            this_header_texts = ""  # start collecting header text

            # assume file is bgzipped
            with bgzip.BGZipReader(raw) as f:
                while True:
                    # read a small block of bytes at a time
                    if not (d := f.read(10 * 1024)):
                        break

                    # concat the latest chunk of text
                    text = d.tobytes().decode()
                    d.release()
                    this_header_texts += text

                    # check if we've reached the end of the header section and get one
                    # more chunk
                    if break_next_time:
                        break
                    elif "\n#CHROM" in this_header_texts:
                        break_next_time = True

            # extract the header lines and the column headers
            this_header_lines = this_header_texts.split("\n")

            if col_header_line is None:
                # extract the line with column names
                col_header_line = [
                    x for x in this_header_lines if x.startswith("#CHROM")
                ][0]

            this_header_lines = [x for x in this_header_lines if x.startswith("##")]

            # add to the collected header lines
            header_lines.extend(this_header_lines)

    # de-dup but keep original order of lines
    return (
        pd.Series([*header_lines, col_header_line])
        .astype("string")
        .drop_duplicates()
        .tolist()
    )


def make_chunks(header_lines: list[str], chunk_size: int) -> TypedDataFrame[Chunk]:
    """
    Collect the contigs in the VCF header lines and their lengths and split them up
    into chunks for batched processing.

    :param header_lines: the list of distinct VCF header lines
    :param chunk_size: number of lines in the VCFs to process at a time
    :return: a data frame of chunked regions
    """

    # contig header lines contain names and lengths of chroms
    contig_lines = [x for x in header_lines if x.startswith("##contig")]
    chrom_lengths_search = [
        re.search(r"^##contig=<ID=chr([^,]+),length=(\d+)>$", x) for x in contig_lines
    ]
    chrom_lengths = pd.DataFrame(
        [x.groups() for x in chrom_lengths_search if x is not None],
        columns=["chrom", "length"],
    ).astype({"chrom": "string", "length": "int64"})

    # remove things like decoy contigs
    valid_chroms = {str(x) for x in range(1, 23)}.union({"X", "Y"})
    chrom_lengths = chrom_lengths.loc[chrom_lengths["chrom"].isin(valid_chroms)]

    # start splitting regions into smaller chunks
    chunks = []

    for _, r in chrom_lengths.iterrows():
        # `bcftools view` expects 1-index regions in TSVs used for splitting
        start = 1

        while start - 1 < r["length"]:
            end = min(start + chunk_size - 1, r["length"].squeeze())
            chunks.append(
                {
                    "chrom": f"chrom{r['chrom']}",
                    "start": start,
                    "end": end,
                    "split_dir_name": f"chrom{r['chrom']}_{start}-{end}",
                }
            )
            start = end + 1

    return type_data_frame(pd.DataFrame(chunks), Chunk)


def chunk_vcfs(
    vcf_paths: list[Path | str], tmp_dir: str, chunks: TypedDataFrame[Chunk]
) -> None:
    """
    Index each VCF file and concurently split them up into chunks using `bcftools`.

    :param vcf_paths: a list of paths to VCF files
    :param tmp_dir: a path to a temporary directory
    :param chunks: a data frame of chunked regions
    """

    # need to make sure each VCF is tabix-indexed first
    with concurrent.futures.ProcessPoolExecutor() as executor:
        logger.info(f"Indexing {len(vcf_paths)} VCFs")
        for path in vcf_paths:
            executor.submit(index_vcf, path)

    # make the output split dirs
    chunks["split_dir"] = chunks["split_dir_name"].apply(
        lambda x: os.path.join(tmp_dir, x)
    )

    for _, r in chunks.iterrows():
        os.makedirs(r["split_dir"].squeeze())

    # cross VCF paths with chunk paths
    vcf_chunks = pd.DataFrame({"path": vcf_paths}).merge(chunks, how="cross")
    vcf_chunks["chunk_path"] = vcf_chunks.apply(
        lambda x: os.path.join(x["split_dir"], os.path.basename(x["path"])), axis=1
    ).str.rstrip(".gz")
    vcf_chunks = vcf_chunks.drop(columns=["split_dir_name", "split_dir"])

    # write the chunks
    with concurrent.futures.ProcessPoolExecutor() as executor:
        logger.info(f"Writing {len(vcf_chunks)} VCF chunks")
        for r in vcf_chunks.to_dict(orient="records"):
            executor.submit(write_vcf_chunk, r)


def index_vcf(path: Path | str) -> None:
    """
    Create a tabix index for a VCF file using `bcftools`.

    :param path: a path to a VCF file
    """

    subprocess.run(["bcftools", "index", path, "--force"])


def write_vcf_chunk(r: dict[Hashable, Any]) -> None:
    """
    Use `bcftools view` to subset a region of a VCF file and output to a file.

    :param r: a dictionary indicating the path and region (chr, start, end) for a chunk
    """

    with open(r["chunk_path"], "w") as f:
        # pipe filtered region of VCF to stdout
        p1 = subprocess.Popen(
            [
                "bcftools",
                "view",
                r["path"],
                "--no-header",
                "--regions",
                f"{r['chrom']}:{r['start']}-{r['end']}",
                "--output-type",
                "v",
            ],
            stdout=subprocess.PIPE,
        )

        # replace invalid UTF characters (e.g. nbsp's from Funcotator)
        p2 = subprocess.Popen(
            ["iconv", "-f", "utf-8", "-t", "utf-8", "-c", "-s"],
            stdin=p1.stdout,
            stdout=f,
        )

        p1.stdout.close()  # pyright: ignore
        _ = p2.communicate(timeout=300)


def process_chunks(
    header_lines: list[str],
    out_path: Path | str,
    tmp_dir: str,
    chunks: TypedDataFrame[Chunk],
    n_vcfs: int,
) -> None:
    """
    Iterate over chunks and output their INFO-merged lines to the output file.

    :param header_lines: the list of distinct VCF header lines
    :param out_path: the path to write the merged VCF file
    :param tmp_dir: a path to a temporary directory
    :param chunks: a data frame of chunked regions
    :param n_vcfs: the number of input VCF files
    """

    with open(out_path, "wb") as raw:
        with bgzip.BGZipWriter(raw) as f:
            # start writing bgzipped output VCF
            logger.info(f"Writing header to {out_path}")
            s = "\n".join(header_lines) + "\n"
            f.write(s.encode())

            for _, r in chunks.iterrows():
                split_dir = os.path.join(tmp_dir, r["split_dir_name"].squeeze())

                # there should be one file per split per input VCF
                split_files = glob(os.path.join(split_dir, "*.vcf"), root_dir=tmp_dir)

                assert len(split_files) == n_vcfs, (
                    f"Expecting {n_vcfs} in {r['split_dir_name']}, "
                    f"but found {len(split_files)}"
                )

                logger.info(
                    f"Merging {len(split_files)} VCF chunks for {r['split_dir_name']}"
                )
                df = merge_chunks(split_files)
                [os.remove(x) for x in split_files]

                logger.info(f"Writing merged {r['split_dir_name']} rows to {out_path}")
                s = df.to_csv(header=False, index=False, sep="\t", quoting=QUOTE_NONE)
                f.write(s.encode())


def merge_chunks(split_files: list[str]) -> TypedDataFrame[Vcf]:
    """
    Reduce a list of VCF files into a single data frame, merging their INFO fields.

    :param split_files: a list of headerless VCF files for this chunk/region
    :return: a data frame of the INFO-merged records
    """

    df = None

    # reduce VCF chunks
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(read_vcf, p) for p in split_files]

        for future in concurrent.futures.as_completed(futures):
            df = merge_two_vcfs(df, future.result())

    assert df is not None

    # de-dup each info field's annotations
    df["info"] = df["info"].apply(collect_uniq_annots).replace({"": "."})

    # df contains only one chrom, so sorting by pos is sufficient
    df = df.sort_values("pos")

    # put back the NA value indicator
    df = df.fillna(".")

    return type_data_frame(df, Vcf)


def read_vcf(path: str) -> TypedDataFrame[Vcf]:
    """
    Read a VCF file as a data frame.

    :param path: a path to a VCF file
    :return: a data frame representation of the VCF's records
    """

    with open(path, "r") as f:
        # noinspection PyTypeChecker
        df = pd.read_csv(
            f,
            sep="\t",
            header=None,
            names=[
                "chrom",
                "pos",
                "id",
                "ref",
                "alt",
                "qual",
                "filter",
                "info",
                "format",
                "values",
            ],
            dtype="string",
            na_values=["."],
            keep_default_na=False,
            quoting=QUOTE_NONE,
            encoding_errors="backslashreplace",
        )

        # remove header rows
        df = df.loc[~df["chrom"].str.startswith("#")]

        return type_data_frame(df, Vcf)


def merge_two_vcfs(
    df1: TypedDataFrame[Vcf] | None, df2: TypedDataFrame[Vcf]
) -> TypedDataFrame[Vcf]:
    """
    Merge two data frames of VCF records on (CHROM, POS, REF, ALT) columns. Columns
    ID, QUAL, FILTER, FORMAT, and VALUES are expected to either by identical in both
    data frames or missing on one side, in which case the value is filled. INFO fields
    are concatenated.

    :param df1: a data frame of VCF records
    :param df2: a data frame of VCF records
    :return: an INFO-merged data frame
    """

    if df1 is None:
        # starting base data frame
        return type_data_frame(df2, Vcf)

    # joining `id` and `info` fields to existing df
    df = df1.merge(
        df2,
        how="outer",
        on=["chrom", "pos", "ref", "alt"],
        suffixes=("", "2"),
    )

    # after correcting for variants potentially only being on one side of the merge,
    # ensure that columns that should be invariant are the same on each side
    for c in ["id", "qual", "filter", "format", "values"]:
        df[c] = df[c].fillna(df[f"{c}2"])
        assert df[c].eq(df[f"{c}2"]).all()

    # concat info fields for each row (will de-dup later)
    df["info"] = df["info"].fillna("") + ";" + df["info2"].fillna("")

    df = df.drop(columns=["id2", "qual2", "filter2", "format2", "values2", "info2"])

    return type_data_frame(df, Vcf)


def collect_uniq_annots(x: str) -> str:
    """
    Dedup a single INFO field's values.

    :param x: an INFO field (semicolon-separated values)
    :return: a deduped INFO field
    """

    annots = {x for x in re.split(r";+", x) if len(x) > 0}
    annots = sorted(list(annots))
    return ";".join(annots)


def type_data_frame(
    df: pd.DataFrame,
    pandera_schema: Type[PanderaBaseSchema],
    remove_unknown_cols: bool = False,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Coerce a data frame into one specified by a Pandera schema and optionally remove
    unknown columns.

    :param df: a data frame
    :param pandera_schema: a Pandera schema
    :param remove_unknown_cols: remove columns not specified in the schema
    :return: a data frame validated with the provided Pandera schema
    """

    if len(df) == 0:
        # make an empty data frame that conforms to the Pandera schema
        s = pandera_schema.to_schema()

        # `example` doesn't know how to instantiate dicts, so do that manually
        dict_cols = []

        for c in s.columns:
            if s.columns[c].dtype.type is dict:
                dict_cols.append(c)
                s = s.remove_columns([c])

        df = pd.DataFrame(s.example(size=0))

        if len(dict_cols) > 0:
            for c in dict_cols:
                df[c] = {}

    elif remove_unknown_cols:
        df_cols = list(pandera_schema.to_schema().columns.keys())
        df = df.loc[:, df_cols]

    return TypedDataFrame[pandera_schema](df)
