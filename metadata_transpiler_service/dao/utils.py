# Copyright 2021 - 2022 Universität Tübingen, DKFZ and EMBL
# for the German Human Genome-Phenome Archive (GHGA)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Utils for DAO"""

import shutil
from tempfile import SpooledTemporaryFile
from typing import IO, Dict, List, Union

import pandas as pd
from pandas import DataFrame

TMP_FILE = "tmp_file.xlsx"


async def save_to_tmp_file(file: Union[SpooledTemporaryFile, IO]):
    """Save File"""
    with open(TMP_FILE, "wb") as buffer:
        shutil.copyfileobj(file, buffer)


async def read_submission_sheets(
    file: Union[SpooledTemporaryFile, IO], sheet_names: List
) -> Dict[str, DataFrame]:
    """Read excel file"""
    await save_to_tmp_file(file)
    sheets = pd.read_excel(
        io=TMP_FILE,
        sheet_name=sheet_names,
        skiprows=list(range(2, 5)),
        header=1,
        index_col=0,
        keep_default_na=False,
        dtype=str,
    )
    return sheets


async def read_value(
    submission_sheet: DataFrame, xls_line_index: int, xls_col_name: str
) -> Union[str, None]:
    """Get a raw value from Excel (as a string)"""
    return submission_sheet.iloc[xls_line_index][xls_col_name]
