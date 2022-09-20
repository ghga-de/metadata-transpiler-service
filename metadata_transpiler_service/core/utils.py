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
Module containing the utils for convertion.
"""

from typing import Dict, List, Union


async def map_to(field: str, mapping_scheme: Dict) -> Union[str, None]:
    """Returns the field name in spreadsheet which contains data entered
    for a field in submission"""
    if field not in mapping_scheme:
        return None
    return mapping_scheme[field]


async def add_unique(list1: List, list2: List) -> List:
    """Add only unique elements to the list

    Args:
        list1 (List): original list
        list2 (List): list to be added

    Returns:
        List: concatenated list
    """
    for elem in list2:
        if elem not in list1:
            list1.append(elem)

    return list1


async def exists_in(alias: str, embedded_list: List) -> bool:
    """Check if an object is already in the list

    Args:
        alias (str): alias of the object to check
        embedded_list (List): the list of objects

    Returns:
        bool: true, if the object is already in the list, otherwise false
    """
    hits = [x for x in embedded_list if x["alias"] == alias]

    if len(hits) != 0:
        return True

    return False


async def get_from_list(alias: str, embedded_list: List) -> Union[Dict, None]:
    """Select an object from embedded list by alias (only the first hit)

    Args:
        alias (str): alias of the object to check
        embedded_list (List): the list of objects

    Returns:
        Dict: the first found object, if not found, then None
    """
    hits = [x for x in embedded_list if x["alias"] == alias]

    if len(hits) == 1:
        return hits[0]

    return None


async def create_list(input_object: object) -> List:
    """Transform object to list, if it is not a list already"""
    if isinstance(input_object, list):
        return input_object
    return [input_object]
