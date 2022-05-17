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

from typing import Dict, Union


async def map_to(field: str, mapping_scheme: Dict) -> Union[str, None]:
    """Returns the field name in spreadsheet which contains data entered
    for a field in submission"""
    if field not in mapping_scheme:
        return None
    return mapping_scheme[field]
