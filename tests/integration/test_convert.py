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

"""Test the api module"""

from zipfile import ZipFile

from fastapi import status
from fastapi.testclient import TestClient

from metadata_transpiler_service.api.main import app

from ..fixtures import BASE_DIR


def test_convert_xlsx():
    """Test the convert endpoint"""

    client = TestClient(app)
    file_path_zip = (
        BASE_DIR / "test_data" / "submission_example" / "submission_example.zip"
    )

    with ZipFile(file_path_zip, "r") as zip_obj:
        xls_file = {"file": zip_obj.open("submission_example.xlsx", "r")}
        response = client.post("/convert", files=xls_file)

        assert response.status_code == status.HTTP_200_OK

        create_submission_entity = response.json()

        assert "has_study" in create_submission_entity
        assert "has_experiment" in create_submission_entity
        assert len(create_submission_entity["has_experiment"]) == 11
        assert "has_project" in create_submission_entity
        assert "has_publication" in create_submission_entity
        assert "has_sample" in create_submission_entity
        assert len(create_submission_entity["has_sample"]) == 22

        create_sample_entity = create_submission_entity["has_sample"]
        assert "case_control_status" in create_sample_entity[0]

        assert "has_individual" in create_submission_entity
        assert "has_protocol" in create_submission_entity
        assert len(create_submission_entity["has_protocol"]) == 22
        assert "has_analysis" in create_submission_entity
        assert "has_file" in create_submission_entity
