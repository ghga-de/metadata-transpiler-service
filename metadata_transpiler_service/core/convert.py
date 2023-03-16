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

"""
Module containing the convertion logic.
"""

from typing import Any, Dict, List, Optional, Union

from pandas import DataFrame

from metadata_transpiler_service.config import CONFIG
from metadata_transpiler_service.core.schema import (
    get_embedded_schema,
    get_schema_by_type,
    is_array,
    is_enum,
    is_integer,
)
from metadata_transpiler_service.core.utils import (
    add_unique,
    create_list,
    exists_in,
    map_to,
)
from metadata_transpiler_service.creation_models import CreateSubmission
from metadata_transpiler_service.dao.utils import read_value

SCHEMA_VERSION = CONFIG.schema_version


async def generate_json_from(submission_sheets: Dict, submission_map: Dict):
    """Generate a submission JSON from the provided spreadsheet (multiple sheets)"""

    submission_object: Dict[str, Any] = {}
    submission_fields: Dict[str, Any]
    sheet_names = list(submission_map.keys())
    submission_schema = await get_schema_by_type("CreateSubmission", None)

    for sheet_name in sheet_names:
        submission_sheet = submission_sheets[sheet_name]
        submission_fields = submission_map[sheet_name]

        for submission_field, embedded_submission_fields in submission_fields.items():
            submission_object = await generate_embedded_json_from(
                submission_sheet,
                submission_field,
                embedded_submission_fields,
                submission_object,
                submission_schema,
            )

    submission_object["schema_type"] = "CreateSubmission"
    submission_object["schema_version"] = SCHEMA_VERSION

    return submission_object


async def generate_embedded_json_from(
    submission_sheet: DataFrame,
    submission_field: str,
    embedded_submission_fields: Dict,
    json_object: Dict,
    submission_schema: Dict,
) -> Dict:
    """Generate JSON for embedded submission objects"""

    if submission_field.startswith("link_"):
        submission_field_name = submission_field.replace("link_", "has_")
    else:
        submission_field_name = submission_field

    if await is_array(submission_schema["properties"][submission_field_name]):
        embedded_list: List = []
        list_to_iterate = []
        ref_dictionary: Dict[str, Dict] = {}
        list_to_iterate = await create_list(embedded_submission_fields)
        for line in range(0, len(submission_sheet.index)):
            for item in list_to_iterate:
                embedded_json_fields = await build_embedded_submission_fields(
                    submission_sheet,
                    item,
                    submission_field_name,
                    line,
                )
                if embedded_json_fields:
                    if submission_field.startswith("link_"):
                        ref_dictionary = await get_refs(
                            embedded_json_fields, ref_dictionary
                        )

                    if not await exists_in(
                        embedded_json_fields["alias"], embedded_list
                    ):
                        embedded_list.append(embedded_json_fields)

        if submission_field.startswith("link_"):
            json_object = await update_object_with_refs(
                json_object, submission_field_name, ref_dictionary
            )
        else:
            json_object[submission_field] = embedded_list

    else:
        json_object[submission_field] = await build_embedded_submission_fields(
            submission_sheet, embedded_submission_fields, submission_field_name
        )
    return json_object


async def build_embedded_submission_fields(
    submission_sheet: Dict, map_embedded_submission_fields: Dict, field: str, count=0
) -> Union[Dict, None]:
    """Build embedded JSON object using data from Spreadsheet"""

    schema = await get_embedded_schema(map_embedded_submission_fields, field)
    field_name = await map_to("alias", map_embedded_submission_fields)
    if not field_name:
        raise ValueError(
            f"No field in the input file maps "
            f"to '{field_name}' -> 'alias' in submission"
        )
    field_value = await read_value(submission_sheet, count, field_name)
    if not field_value:
        return None

    embedded_submission_fields: Dict = {}
    for json_field, xls_col_name in map_embedded_submission_fields.items():
        if json_field.startswith("_"):
            continue
        json_field = json_field.replace("_additional_attribute", "_attribute")
        if isinstance(xls_col_name, list):
            if (json_field not in embedded_submission_fields) or (
                not embedded_submission_fields[json_field]
            ):
                embedded_submission_fields[json_field] = []
            for single_col_name in xls_col_name:
                field_value = await read_value(submission_sheet, count, single_col_name)
                final_value = await convert_value(
                    field_value, json_field, schema, single_col_name
                )
                if final_value:
                    embedded_submission_fields[json_field] = (
                        embedded_submission_fields[json_field] + final_value
                    )
        else:
            field_value = await read_value(submission_sheet, count, xls_col_name)
            if (json_field not in embedded_submission_fields) or (
                not embedded_submission_fields[json_field]
            ):
                embedded_submission_fields[json_field] = await convert_value(
                    field_value, json_field, schema
                )
            else:
                embedded_submission_fields[json_field] = embedded_submission_fields[
                    json_field
                ] + await convert_value(field_value, json_field, schema)

    embedded_submission_fields["schema_type"] = schema["title"]
    embedded_submission_fields["schema_version"] = SCHEMA_VERSION

    return embedded_submission_fields


async def get_refs(embedded_submission_fields: Dict, ref_dictionary: Dict) -> Dict:
    """Add references from an object to the reference dictionary"""
    parent = embedded_submission_fields["alias"]
    for field_name in embedded_submission_fields:
        if field_name.startswith("has_"):
            ref = embedded_submission_fields[field_name]
            if field_name not in ref_dictionary:
                ref_dictionary[field_name] = {}
            if parent in ref_dictionary[field_name]:
                ref_dictionary[field_name][parent] = await add_unique(
                    ref_dictionary[field_name][parent], ref
                )
            else:
                ref_dictionary[field_name][parent] = ref

    return ref_dictionary


async def update_object_with_refs(
    submission_object: Dict, field_name: str, ref_dictionary: Dict
) -> Dict:
    """Add references to submission object for a field"""
    ref_list = []
    for obj in submission_object[field_name]:
        for link_key, linked_obj in ref_dictionary.items():
            obj[link_key] = linked_obj[obj["alias"]]
        ref_list.append(obj)
    submission_object[field_name] = ref_list
    return submission_object


async def convert_value(
    field_value: Optional[str],
    field_name: str,
    model_schema: Dict,
    single_col_name: Optional[str] = None,
) -> Union[Dict, List, str, int, None]:
    """Normalize raw field value and transform to requested format"""
    if not field_value:
        return None
    if await is_enum(model_schema["properties"][field_name]):
        field_value = await normalize(field_value)
    if await is_integer(model_schema["properties"][field_name]):
        return int(float(field_value))
    final_value = await transform(
        field_value, field_name, model_schema["title"], single_col_name
    )
    if not await is_array(model_schema["properties"][field_name]):
        return final_value[0]

    return final_value


async def normalize(cell_value: str) -> str:
    """Normalize value"""
    return cell_value.replace(" ", "_").lower()


async def transform(
    cell_value: str, cell_type: str, parent: str, single_col_name: Optional[str] = None
) -> List:
    """Transform string value to list of values depending on the required type
    (string, array, attribute(s), enumeration, embedded entity (concept))"""
    list_values = cell_value.split(";")
    if cell_type == "has_attribute":
        list_attr = []
        if single_col_name is not None:
            obj_attr = {"key": single_col_name, "value": cell_value}
            list_attr.append(obj_attr)
        else:
            for attr in list_values:
                pair = attr.split("=")
                obj_attr = {
                    "key": pair[0],
                    "value": pair[1],
                }
                list_attr.append(obj_attr)
        return list_attr
    if cell_type.startswith("has_") and cell_type not in CreateSubmission.__fields__:
        schema = await get_schema_by_type(parent, cell_type)
        if "concept_name" in schema["properties"]:
            list_entities = []
            for entity in list_values:
                obj_entity = {
                    "concept_name": entity,
                    "schema_type": schema["title"],
                    "schema_version": SCHEMA_VERSION,
                }
                list_entities.append(obj_entity)
            return list_entities
    return list_values
