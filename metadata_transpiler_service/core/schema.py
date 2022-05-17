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
Module containing the utils for the model scheme.
"""

from importlib import import_module
from typing import Dict, Union

MODELS_MODULE_NAME = "metadata_transpiler_service.creation_models"


async def is_array(embedded_submission_fieldsect: Dict):
    """Returns true if an embedded object is an array"""
    if "type" in embedded_submission_fieldsect:
        if embedded_submission_fieldsect["type"] == "array":
            return True
    if "anyOf" in embedded_submission_fieldsect:
        if await is_array(embedded_submission_fieldsect["anyOf"][0]):
            return True
    return False


async def is_enum(embedded_submission_fieldsect: Dict):
    """Returns true if an embedded object is of type enumeration"""
    if "allOf" in embedded_submission_fieldsect and embedded_submission_fieldsect[
        "allOf"
    ][0]["$ref"].endswith("Enum"):
        return True
    return False


async def get_reference(schema: Dict) -> Dict:
    """Extract reference to scheme if exists, otherwise the original schema"""
    if "$ref" in schema:
        class_name = schema["$ref"].split("/definitions/")[1]
        schema = schema["definitions"][class_name]
    return schema


async def get_schema(schema_type: str, field_name: Union[str, None]) -> Dict:
    """
    Get schema for a metadata object. If field is provided, get embedded schema for a field
    """
    module = import_module(MODELS_MODULE_NAME)
    new_class = getattr(module, schema_type)
    main_schema = new_class.schema()
    if not field_name:
        return main_schema

    embedded_schema = await get_reference(
        new_class.__fields__[field_name].sub_fields[0].type_.schema()
    )
    return embedded_schema
