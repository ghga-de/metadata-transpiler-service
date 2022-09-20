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
Module containing the inference methods.
"""

from typing import Dict, List

from metadata_transpiler_service.core.schema import get_schema_by_type, is_array
from metadata_transpiler_service.core.utils import (
    add_unique,
    create_list,
    get_from_list,
)


async def infer_missing_fields(submission: Dict) -> Dict:
    """Infer missing field in the submission JSON"""
    submission = await add_alias("study", "dataset", submission)
    submission = await add_alias("study", "experiment", submission)
    submission = await add_alias("study", "analysis", submission)
    submission = await autofill_alias("sample", "dataset", "experiment", submission)
    submission = await autofill_alias("file", "dataset", "experiment", submission)

    return submission


async def add_alias(embedded_object: str, parent_object: str, submission: Dict) -> Dict:
    """Add missing alias to embedded object from parent object"""

    field_name = "has_" + parent_object
    schema_name = "Create" + parent_object.capitalize()
    model_schema = await get_schema_by_type(schema_name, None)

    embedded_field_name = "has_" + embedded_object
    if embedded_field_name not in submission:
        return submission
    alias = submission[embedded_field_name]["alias"]
    if field_name not in submission:
        return submission
    object_to_change = submission[field_name]
    if isinstance(object_to_change, list):
        object_list = []
        for single_object in object_to_change:
            if await is_array(model_schema["properties"][embedded_field_name]):
                single_object[embedded_field_name] = [alias]
            else:
                single_object[embedded_field_name] = alias
            object_list.append(single_object)
        submission[field_name] = object_list
    else:
        if await is_array(model_schema["properties"][embedded_field_name]):
            object_to_change[embedded_field_name] = [alias]
        else:
            object_to_change[embedded_field_name] = alias

        submission[field_name] = object_to_change

    return submission


async def autofill_alias(
    embedded_object: str, parent_object: str, source_object: str, submission: Dict
) -> Dict:
    """Collect list of aliases for embedded object from source object
    to autofill the missing entries in parent object"""

    field_name = "has_" + parent_object
    embedded_field_name = "has_" + embedded_object
    source_field_name = "has_" + source_object

    if source_field_name not in submission:
        return submission

    list_to_iterate = await create_list(submission[field_name])
    object_list = []
    for single_object in list_to_iterate:
        if (
            single_object[embedded_field_name] is None
            or len(single_object[embedded_field_name]) == 0
        ):
            if source_field_name in single_object:
                list_of_source_ids = single_object[source_field_name]
                list_of_embedded_ids: List[str] = []
                for item in list_of_source_ids:
                    extract_from = await get_from_list(
                        item, submission[source_field_name]
                    )
                    if extract_from is None or embedded_field_name not in extract_from:
                        continue
                    list_of_embedded_ids = await add_unique(
                        list_of_embedded_ids, extract_from[embedded_field_name]
                    )
            single_object[embedded_field_name] = list_of_embedded_ids
        object_list.append(single_object)
    if isinstance(submission[field_name], list):
        submission[field_name] = object_list
    else:
        submission[field_name] = object_list[0]

    return submission
