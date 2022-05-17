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
Mapping for transpiler.
"""

from typing import Any, Dict

submission_map: Dict[str, Any] = {
    "Core_Submission_Sheet": {
        "has_project": {
            "alias": "submitter_project_alias",
            "title": "project_title",
            "description": "project_description",
            "has_attribute": "project_attributes",
        },
        "has_study": {
            "alias": "submitter_study_alias",
            "title": "study_title",
            "description": "study_description",
            "type": "study_type",
            "affiliation": "study_affiliations",
            "has_attribute": "study_attributes",
            "has_project": "submitter_project_alias",
            "has_publication": "publication_id",
        },
        "has_experiment": {
            "alias": "submitter_experiment_alias",
            "title": "experiment_title",
            "description": "experiment_description",
            "biological_replicates": "experiment_biological_replicates",
            "technical_replicates": "experiment_technical_replicates",
            "experimental_replicates": "experiment_experimental_replicates",
        },
        "has_publication": {
            "alias": "publication_id",
            "title": "publication_title",
            "abstract": "publication_abstract",
            "xref": "publication_xref",
        },
    },
    "Sample_Submission_Sheet": {
        "has_sample": {
            "alias": "submitter_sample_alias",
            "name": "sample_name",
            "description": "sample_description",
            "type": "sample_type",
            "vital_status_at_sampling": "sample_vital_status_at_sampling",
            "isolation": "sample_isolation",
            "storage": "sample_storage",
            "xref": "sample_xref",
            "has_attribute": "sample_attributes",
            "has_anatomical_entity": "sample_tissue",
            "has_individual": "submitter_individual_alias",
            "has_biospecimen": "submitter_biospecimen_alias",
        },
        "link_experiment": {
            "alias": "submitter_experiment_alias",
            "has_sample": "submitter_sample_alias",
        },
        "has_biospecimen": {
            "alias": "submitter_biospecimen_alias",
            "xref": "biospecimen_external_id",
            "name": "biospecimen_name",
            "description": "biospecimen_description",
            "isolation": "biospecimen_isolation",
            "storage": "biospecimen_storage",
            "has_individual": "submitter_individual_alias",
        },
        "has_individual": {
            "alias": "submitter_individual_alias",
            "sex": "individual_sex",
            "age": "individual_age",
            "has_phenotypic_feature": "individual_phenotypic_feature",
            "vital_status": "individual_vital_status",
            "karyotype": "individual_karyotype",
            "geographical_region": "individual_geographical_region",
        },
    },
    "Experimental_Submission_Sheet": {
        "has_protocol": [
            {
                "_model": "CreateLibraryPreparationProtocol",
                "alias": "submitter_libprep_alias",
                "description": "libprep_description",
                "library_name": "libprep_library_name",
                "library_layout": "libprep_library_layout",
                "library_type": "libprep_library_type",
                "library_selection": "libprep_library_selection",
                "library_preparation": "libprep_library_preparation",
                "library_preparation_kit_retail_name": "libprep_kit_retail_name",
                "library_preparation_kit_manufacturer": "libprep_kit_manufacturer",
                "target_regions": "libprep_target_regions",
                "rnaseq_strandedness": "libprep_rnaseq_strandedness",
                "primer": "libprep_primer",
                "end_bias": "libprep_end_bias",
                "has_attribute": "libprep_attributes",
            },
            {
                "_model": "CreateSequencingProtocol",
                "alias": "submitter_sequencing_alias",
                "description": "seq_description",
                "type": "seq_type",
                "instrument_model": "seq_instrument_model",
                "sequencing_center": "seq_center",
                "index_sequence": "seq_index_sequence",
                "sequencing_read_length": "seq_read_length",
                "paired_or_single_end": "seq_paired_or_single_end",
                "target_coverage": "seq_target_coverage",
                "lane_number": "seq_lane_number",
                "flow_cell_id": "seq_flow_cell_id",
                "flow_cell_type": "seq_flow_cell_type",
                "umi_barcode_read": "seq_umi_barcode_read",
                "umi_barcode_offset": "seq_umi_barcode_offset",
                "umi_barcode_size": "seq_umi_barcode_size",
                "cell_barcode_read": "seq_cell_barcode_read",
                "cell_barcode_offset": "seq_cell_barcode_offset",
                "cell_barcode_size": "seq_cell_barcode_size",
                "sample_barcode_read": "seq_sample_barcode_read",
                "has_attribute": "seq_attributes",
            },
        ],
        "link_experiment": {
            "alias": "submitter_experiment_alias",
            "has_protocol": [
                "submitter_libprep_alias",
                "submitter_sequencing_alias",
            ],
            "has_file": "submitter_file_alias",
        },
    },
    "Analysis_Submission_Sheet": {
        "has_analysis": {
            "alias": "submitter_analysis_alias",
            "description": "analysis_description",
            "type": "analysis_type",
            "reference_genome": "analysis_reference_genome",
            "reference_chromosome": "analysis_reference_chromosome",
            "has_study": "submitter_study_alias",
        },
        "link_analysis": {
            "alias": "submitter_analysis_alias",
            "has_input": "submitter_input_file_alias",
            "has_output": "submitter_output_file_alias",
        },
    },
    "File_Submission_Sheet": {
        "has_file": {
            "alias": "submitter_file_alias",
            "name": "file_name",
            "checksum": "file_checksum",
            "checksum_type": "file_checksum_type",
            "format": "file_format",
        }
    },
}
