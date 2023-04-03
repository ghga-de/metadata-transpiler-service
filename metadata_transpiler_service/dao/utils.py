# Copyright 2021 - 2023 Universität Tübingen, DKFZ, EMBL, and Universität zu Köln
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

import json
import os
import shutil
from tempfile import SpooledTemporaryFile
from typing import IO, Dict, List, Union

import pandas as pd
from fastapi import HTTPException
from pandas import DataFrame

from metadata_transpiler_service import BASE_DIR
from metadata_transpiler_service.config import CONFIG

MAPPING_URL = CONFIG.mapping_url
TMP_FILE = "tmp_file.xlsx"

SKIP_ROWS = list(range(2, 6))
HEADER = 1
INDEX_COL = 0


async def save_to_tmp_file(file: Union[SpooledTemporaryFile, IO]):
    """Save File"""
    with open(TMP_FILE, "wb") as buffer:
        shutil.copyfileobj(file, buffer)


async def delete_tmp_file():
    """Delete temporal file"""
    if os.path.exists(TMP_FILE):
        os.remove(TMP_FILE)


async def read_mapping_file() -> Dict:
    """Read mapping JSON from file"""
    file_path = str(BASE_DIR) + MAPPING_URL
    try:
        with open(file_path, "r", encoding="utf8") as file:
            submission_map = json.load(file)
    except FileNotFoundError as fnfe:
        raise HTTPException(
            status_code=404,
            detail=(f"Cannot find the mapping file '{file_path}': " f"{fnfe}"),
        ) from fnfe
    return submission_map


async def read_submission_sheets(
    file: Union[SpooledTemporaryFile, IO], sheet_names: List
) -> Dict[str, DataFrame]:
    """Read excel file"""
    await save_to_tmp_file(file)
    sheets = pd.read_excel(
        io=TMP_FILE,
        sheet_name=sheet_names,
        skiprows=SKIP_ROWS,
        header=HEADER,
        index_col=INDEX_COL,
        keep_default_na=False,
        dtype=str,
    )
    await delete_tmp_file()
    return sheets


async def read_value(
    submission_sheet: DataFrame, xls_line_index: int, xls_col_name: str
) -> Union[str, None]:
    """Get a raw value from Excel (as a string)"""
    return submission_sheet.iloc[xls_line_index][xls_col_name]
