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

"""
Module containing the main FastAPI router and (optionally) top-level API enpoints.
Additional endpoints might be structured in dedicated modules
(each of them having a sub-router).
"""

from tempfile import SpooledTemporaryFile
from typing import IO, Union

from fastapi import FastAPI, File, HTTPException, UploadFile
from ghga_service_chassis_lib.api import configure_app

from metadata_transpiler_service.config import CONFIG
from metadata_transpiler_service.core.convert import generate_json_from
from metadata_transpiler_service.core.infer import infer_missing_fields
from metadata_transpiler_service.creation_models import CreateSubmission
from metadata_transpiler_service.dao.utils import (
    read_mapping_file,
    read_submission_sheets,
)

app = FastAPI()
configure_app(app, config=CONFIG)


@app.get("/", summary="")
async def index():
    """Index"""
    return "Index of the GHGA Metadata Transpiler Service"


@app.post(
    "/convert",
    summary="Given XLSX file, converts it in JSON in Submission format",
    response_model=CreateSubmission,
)
async def convert_xlsx_to_json(file: UploadFile = File(...)):
    """Convert the uploaded spreadsheet into JSON according to the CreateSubmission model"""

    submission_map = await read_mapping_file()
    sheet_names = list(submission_map.keys())
    temp_file: Union[SpooledTemporaryFile, IO] = file.file

    try:
        submission_sheets = await read_submission_sheets(temp_file, sheet_names)
    except Exception as exp:
        raise HTTPException(
            status_code=404,
            detail=(f"Cannot read the input file '{file.filename}': " f"{exp}"),
        ) from exp

    submission_json = await generate_json_from(submission_sheets, submission_map)
    submission_json = await infer_missing_fields(submission_json)

    return submission_json
