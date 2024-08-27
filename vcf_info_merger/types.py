from typing import TypeVar

import pandas as pd
import pandera as pa
from pandera.typing import DataFrame as PaDataFrame
from pandera.typing import Series


class Vcf(pa.DataFrameModel):
    chrom: Series[pd.StringDtype]
    pos: Series[pd.Int64Dtype]
    id: Series[pd.StringDtype] = pa.Field(nullable=True)
    ref: Series[pd.StringDtype] = pa.Field(nullable=True)
    alt: Series[pd.StringDtype] = pa.Field(nullable=True)
    qual: Series[pd.StringDtype] = pa.Field(nullable=True)
    filter: Series[pd.StringDtype] = pa.Field(nullable=True)
    info: Series[pd.StringDtype] = pa.Field(nullable=True)
    format: Series[pd.StringDtype] = pa.Field(nullable=True)
    values: Series[pd.StringDtype] = pa.Field(nullable=True)

    class Config:  # pyright: ignore
        coerce = True
        unique = ["chrom", "pos", "id", "ref", "alt"]


class Chunk(pa.DataFrameModel):
    chrom: Series[pd.StringDtype]
    start: Series[pd.Int64Dtype]
    end: Series[pd.Int64Dtype]
    split_dir_name: Series[pd.StringDtype] = pa.Field(unique=True)

    class Config:  # pyright: ignore
        coerce = True
        unique = ["chrom", "start", "end"]

    @classmethod
    @pa.dataframe_check
    def not_empty(cls, df: pd.DataFrame) -> bool:
        return len(df) > 0


TypedDataFrame = PaDataFrame
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=pa.DataFrameModel)
