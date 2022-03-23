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

"""Defines all dataclasses/classes pertaining to a data model or schema"""


# pylint: disable=too-many-lines

from __future__ import annotations

from enum import Enum
from typing import List, Optional, Union

from pydantic import BaseModel, Field


class CaseControlEnum(str, Enum):
    """
    Enum to capture whether a Sample is to be considered as Case or Control.
    """

    CONTROL = "control"
    CASE = "case"


class BiologicalSexEnum(str, Enum):
    """
    The biological sex of an Individual as determined by their chromosomes.
    """

    XX = "XX"
    XY = "XY"
    NONE = "none"


class UserRoleEnum(str, Enum):
    """
    Enum to capture the different roles a GHGA User can have.
    """

    DATA_REQUESTER = "data requester"
    DATA_STEWARD = "data steward"


class VitalStatusEnum(str, Enum):
    """
    Enum to capture the vital status of an individual.
    """

    ALIVE = "alive"
    DECEASED = "deceased"
    UNKNOWN = "unknown"


class StudyTypeEnum(str, Enum):
    """
    Enum to capture the type of a study.
    """

    WHOLE_GENOME_SEQUENCING = "whole_genome_sequencing"
    METAGENOMICS = "metagenomics"
    TRANSCRIPTOME_ANALYSIS = "transcriptome_analysis"
    RESEQUENCING = "resequencing"
    EPIGENETICS = "epigenetics"
    SYNTHETIC_GENOMICS = "synthetic_genomics"
    FORENSIC_PALEO_GENOMICS = "forensic_paleo_genomics"
    GENE_REGULATION = "gene_regulation"
    CANCER_GENOMICS = "cancer_genomics"
    POPULATION_GENOMICS = "population_genomics"
    RNA_SEQ = "rna_seq"
    EXOME_SEQUENCING = "exome_sequencing"
    POOLED_CLONE_SEQUENCING = "pooled_clone_sequencing"
    OTHER = "other"


class FileTypeEnum(str, Enum):
    """
    Enum to capture file types.
    """

    BAM = "bam"
    COMPLETE_GENOMICS = "complete_genomics"
    CRAM = "cram"
    FASTA = "fasta"
    FASTQ = "fastq"
    PACBIO_HDF5 = "pacbio_hdf5"
    SFF = "sff"
    SRF = "srf"
    VCF = "vcf"


class SubmissionStatusEnum(str, Enum):
    """
    Enum to capture the status of a Submission.
    """

    IN_PROGRESS = "in progress"
    COMPLETED = "completed"


class ReleaseStatusEnum(str, Enum):
    """
    Enum to capture the release status of an entity.
    """

    UNRELEASED = "unreleased"
    RELEASED = "released"


class Attribute(BaseModel):
    """
    A key/value pair that further characterizes an entity.
    """

    key: str = Field(None, description="The key for an attribute.")
    key_type: Optional[str] = Field(
        None,
        description="A semantic type that characterizes the attribute key. \
            Usually this is a term from an ontology. For example, 'MAXO:0000616' \
            indicates that the attribute is a measurement of oxygen saturation in \
            the blood.",
    )
    value: str = Field(
        None,
        description="The value for an attribute. Usually this is a numerical value \
            (without the units).",
    )
    value_type: Optional[str] = Field(
        None,
        description="The value type that characterizes the attribute value. Usually \
            this is a term from an ontology that describes how to interpret the \
            value. For example, 'SIO:001413' indicates that the value is to be \
            interpreted as a percentage.",
    )


class CreateProject(BaseModel):
    """
    Any specifically defined piece of work that is undertaken or attempted to meet
    a single requirement.
    """

    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="Custom attributes for the Project  (eg: Cancer - Colon cancer, \
            prostrate cancer, blood cancer etc)",
    )
    title: str = Field(None, description="Comprehensive title for the project.")
    description: str = Field(
        None,
        description="Short textual description of the project   (Some information on \
            the protocol, sample used and collected etc)  ",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Project(BaseModel):
    """
    Any specifically defined piece of work that is undertaken or attempted to meet
    a single requirement.
    """

    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="Custom attributes for the Project  (eg: Cancer - Colon cancer, \
            prostrate cancer, blood cancer etc)",
    )
    title: str = Field(None, description="Comprehensive title for the project.")
    description: str = Field(
        None,
        description="Short textual description of the project   (Some information on \
            the protocol, sample used and collected etc)  ",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateStudy(BaseModel):
    """
    Studies are experimental investigations of a particular phenomenon. It involves
    a detailed examination and analysis of a subject to learn more about the
    phenomenon being studied.
    """

    affiliation: List[str] = Field(
        None, description="Institutions associated with this study."
    )
    has_experiment: Optional[Union[List[str], List[CreateExperiment]]] = Field(
        None, description="One or more Experiment entities associated with this Study."
    )
    has_analysis: Optional[Union[List[str], List[CreateAnalysis]]] = Field(
        None, description="One or more Analysis entities associated with this Study."
    )
    has_project: Optional[Union[str, CreateProject]] = Field(
        None, description="The project associated with this Study."
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    has_publication: Optional[Union[List[str], List[CreatePublication]]] = Field(
        None, description="One or more Publication entities associated with this Study."
    )
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="Custom key/value pairs that further characterizes the Study. \
            (e.g.: approaches - single-cell, bulk etc)",
    )
    title: str = Field(None, description="Comprehensive title for the study.")
    description: str = Field(
        None,
        description="A detailed description (abstract) that describes the goals of \
            this Study.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    type: Optional[str] = Field(None, description="The type of the Study.")


class Study(BaseModel):
    """
    Studies are experimental investigations of a particular phenomenon. It involves
    a detailed examination and analysis of a subject to learn more about the
    phenomenon being studied.
    """

    affiliation: List[str] = Field(
        None, description="Institutions associated with this study."
    )
    has_experiment: Optional[Union[List[str], List[Experiment]]] = Field(
        None, description="One or more Experiment entities associated with this Study."
    )
    has_analysis: Optional[Union[List[str], List[Analysis]]] = Field(
        None, description="One or more Analysis entities associated with this Study."
    )
    has_project: Optional[Union[str, Project]] = Field(
        None, description="The project associated with this Study."
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    has_publication: Optional[Union[List[str], List[Publication]]] = Field(
        None, description="One or more Publication entities associated with this Study."
    )
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="Custom key/value pairs that further characterizes the Study. \
            (e.g.: approaches - single-cell, bulk etc)",
    )
    status: Optional[ReleaseStatusEnum] = Field(
        None,
        description="The status of a Study. For example, 'released' or 'unreleased'.",
    )
    release_date: Optional[str] = Field(
        None,
        description="The timestamp (in ISO 8601 format) when the entity was released \
            for public consumption.",
    )
    title: str = Field(None, description="Comprehensive title for the study.")
    description: str = Field(
        None,
        description="A detailed description (abstract) that describes the goals of \
            this Study.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    type: Optional[str] = Field(None, description="The type of the Study.")
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateExperiment(BaseModel):
    """
    An experiment is an investigation that consists of a coordinated set of actions
    and observations designed to generate data with the goal of verifying, falsifying,
    or establishing the validity of a hypothesis.
    """

    biological_replicates: Optional[str] = Field(
        None,
        description="A biological replicate is a replicate role that consists of \
            independent biological replicates made from different \
            individual biosamples.",
    )
    technical_replicates: Optional[str] = Field(
        None,
        description="A technical replicate is a replicate role where the same \
            BioSample is use e.g.  the same pool of RNA used to assess technical \
            (as opposed to biological)  variation within an experiment.",
    )
    experimental_replicates: Optional[str] = Field(
        None,
        description="The replicate number of the assay, i.e. the numeric iteration \
            for the assay that was repeated.",
    )
    has_study: Union[str, CreateStudy] = Field(
        None, description="The Study entity associated with this Experiment."
    )
    has_sample: Union[str, CreateSample] = Field(
        None, description="The Sample entity associated with this Experiment."
    )
    has_technology: Optional[Union[str, CreateTechnology]] = Field(
        None, description="The Technology entity associated with this Experiment."
    )
    has_file: Optional[Union[List[str], List[CreateFile]]] = Field(
        None,
        description="One or more Files entities that are generated as output of \
            this Experiment.",
    )
    has_experiment_process: Optional[
        Union[List[str], List[CreateExperimentProcess]]
    ] = Field(
        None,
        description="One or more Experiment Processes entities associated with \
            this Experiment.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    title: str = Field(
        None, description="Name for the experiment (eg: GHGAE_PBMC_RNAseq)."
    )
    description: str = Field(
        None, description="A detailed description of the Experiment."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    type: Optional[str] = Field(None, description="The type of Experiment.")


class Experiment(BaseModel):
    """
    An experiment is an investigation that consists of a coordinated set of actions
    and observations designed to generate data with the goal of verifying, falsifying,
    or establishing the validity of a hypothesis.
    """

    biological_replicates: Optional[str] = Field(
        None,
        description="A biological replicate is a replicate role that consists of \
            independent biological replicates made from different \
            individual biosamples.",
    )
    technical_replicates: Optional[str] = Field(
        None,
        description="A technical replicate is a replicate role where the same \
            BioSample is use e.g.  the same pool of RNA used to assess technical \
            (as opposed to biological)  variation within an experiment.",
    )
    experimental_replicates: Optional[str] = Field(
        None,
        description="The replicate number of the assay, i.e. the numeric iteration \
            for the assay that was repeated.",
    )
    has_study: Union[str, Study] = Field(
        None, description="The Study entity associated with this Experiment."
    )
    has_sample: Union[str, Sample] = Field(
        None, description="The Sample entity associated with this Experiment."
    )
    has_technology: Optional[Union[str, Technology]] = Field(
        None, description="The Technology entity associated with this Experiment."
    )
    has_file: Optional[Union[List[str], List[File]]] = Field(
        None,
        description="One or more Files entities that are generated as output of \
            this Experiment.",
    )
    has_experiment_process: Optional[Union[List[str], List[ExperimentProcess]]] = Field(
        None,
        description="One or more Experiment Processes entities associated with \
            this Experiment.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    title: str = Field(
        None, description="Name for the experiment (eg: GHGAE_PBMC_RNAseq)."
    )
    description: str = Field(
        None, description="A detailed description of the Experiment."
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    type: Optional[str] = Field(None, description="The type of Experiment.")
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateExperimentProcess(BaseModel):
    """
    An Experiment Process is a process that describes how a Sample is transformed to
    a File via an assay. The Experiment Process also keeps track of the Protocol
    used and the Agent that is running the experiment.
    """

    title: Optional[str] = Field(
        None,
        description="A descriptive title that explains the step(s) involved in \
            performing the experiment leading up to the sequencing of the sample \
            and generation of raw data from the instrument. \
            (eg: Sample extraction -> Target Enrichment) ",
    )
    has_input: Optional[Union[str, CreateSample]] = Field(
        None,
        description="The input to the Experiment Process. Usually a Sample entity.",
    )
    has_protocol: Optional[Union[str, CreateProtocol]] = Field(
        None, description="The Protocol entity used by this Experiment Process."
    )
    has_agent: Optional[Union[str, CreateAgent]] = Field(
        None,
        description="The Agent - a software, institution, or human - that is \
            executing or responsible for executing the Experiment Process.",
    )
    has_output: Optional[Union[str, CreateFile]] = Field(
        None,
        description="The output of this Experiment Process. Usually a File entity.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class ExperimentProcess(BaseModel):
    """
    An Experiment Process is a process that describes how a Sample is transformed to
    a File via an assay. The Experiment Process also keeps track of the Protocol
    used and the Agent that is running the experiment.
    """

    title: Optional[str] = Field(
        None,
        description="A descriptive title that explains the step(s) involved in \
            performing the experiment leading up to the sequencing of the sample \
            and generation of raw data from the instrument. \
            (eg: Sample extraction -> Target Enrichment) ",
    )
    has_input: Optional[Union[str, Sample]] = Field(
        None,
        description="The input to the Experiment Process. Usually a Sample entity.",
    )
    has_protocol: Optional[Union[str, Protocol]] = Field(
        None, description="The Protocol entity used by this Experiment Process."
    )
    has_agent: Optional[Union[str, Agent]] = Field(
        None,
        description="The Agent - a software, institution, or human - that is \
            executing or responsible for executing the Experiment Process.",
    )
    has_output: Optional[Union[str, File]] = Field(
        None,
        description="The output of this Experiment Process. Usually a File entity.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateProtocol(BaseModel):
    """
    A plan specification which has sufficient level of detail and quantitative
    information to communicate it between investigation agents, so that different
    investigation agents will reliably be able to independently reproduce the process.
    """

    name: Optional[str] = Field(
        None,
        description="Name of the protocol (eg: Sample extraction_PCR amplification).",
    )
    description: Optional[str] = Field(
        None, description="Detailed description of the Protocol."
    )
    url: Optional[str] = Field(
        None, description="URL for the resource that describes this Protocol."
    )
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="One or more attributes that further characterizes this Protocol.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None,
        description="One or more cross-references for this Protocol. \
            (Eg: manufacturer protocol, protocol from publication etc )",
    )


class Protocol(BaseModel):
    """
    A plan specification which has sufficient level of detail and quantitative
    information to communicate it between investigation agents, so that different
    investigation agents will reliably be able to independently reproduce the process.
    """

    name: Optional[str] = Field(
        None,
        description="Name of the protocol (eg: Sample extraction_PCR amplification).",
    )
    description: Optional[str] = Field(
        None, description="Detailed description of the Protocol."
    )
    url: Optional[str] = Field(
        None, description="URL for the resource that describes this Protocol."
    )
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="One or more attributes that further characterizes this Protocol.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None,
        description="One or more cross-references for this Protocol. \
            (Eg: manufacturer protocol, protocol from publication etc )",
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateLibraryPreparationProtocol(Protocol):
    """
    Information about the library preparation of an Experiment.
    """

    library_name: str = Field(
        None,
        description="A short name identifying the library to potential users. \
            The same name may refer to multiple versions of the same continually \
            updated library.",
    )
    library_layout: str = Field(
        None,
        description="Describe whether the library was sequenced in single-end \
            (forward or reverse) or paired-end mode",
    )
    library_type: str = Field(
        None,
        description="Describe the level of omics analysis (eg: Metagenome, \
            transcriptome, etc) ",
    )
    library_selection: str = Field(
        None,
        description="Whether any method was used to select for or against, enrich, \
            or screen the material being sequenced. Library Selection method \
            (e.g. random, PCA, cDNA, etc )",
    )
    library_construction: str = Field(
        None,
        description="The name of a library construction approach being used. \
            (eg: '10X v2 sequencing' or 'Smart-seq2')",
    )
    library_preparation: str = Field(
        None,
        description="The general method for sequencing library construction \
            (e.g. KAPA PCR-free).",
    )
    library_level: Optional[str] = Field(
        None, description="Single Cell Sequencing or Bulk Sequencing"
    )
    library_construction_kit_retail_name: str = Field(
        None,
        description="A unique identifier for the kit used to construct a \
            genomic library.  This may include the vendor name, kit name and \
            kit version  (e.g. Agilent sure select Human Exome V8, \
            Twist RefSeq Exome, etc.)",
    )
    library_construction_kit_manufacturer: str = Field(
        None, description="Manufacturer of library construction kit"
    )
    primer: Optional[str] = Field(
        None,
        description="The type of primer used for reverse transcription, e.g. \
        'oligo-dT' or 'random' primer. This allows users to identify content of \
        the cDNA library input e.g. enriched for mRNA.",
    )
    end_bias: Optional[str] = Field(
        None,
        description="The end of the cDNA molecule that is preferentially sequenced, \
            e.g. 3/5 prime tag or end, or the full-length transcript.",
    )
    target_regions: Optional[str] = Field(
        None,
        description="Subset of genes or specific regions of the genome, which are \
            most likely to be involved in the phenotype under study. ",
    )
    rnaseq_strandedness: Optional[str] = Field(
        None,
        description="The strandedness of the library, whether reads come from both \
            strands of the cDNA  or only from the first (antisense) or the second \
            (sense) strand.",
    )
    name: str = Field(
        None,
        description="Name of the library preparation protocol (eg: mRNA-seq \
            library preparation).",
    )
    description: str = Field(
        None,
        description="Description about how a sequencing library was prepared \
            (eg: Library construction method).",
    )
    url: Optional[str] = Field(None, description="A URL to a resource.")
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="One or more attributes that further characterizes this \
            Library Preparation Protocol.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class LibraryPreparationProtocol(Protocol):
    """
    Information about the library preparation of an Experiment.
    """

    library_name: str = Field(
        None,
        description="A short name identifying the library to potential users. \
            The same name may refer to multiple versions of the same continually \
            updated library.",
    )
    library_layout: str = Field(
        None,
        description="Describe whether the library was sequenced in single-end \
            (forward or reverse) or paired-end mode",
    )
    library_type: str = Field(
        None,
        description="Describe the level of omics analysis (eg: Metagenome, \
            transcriptome, etc) ",
    )
    library_selection: str = Field(
        None,
        description="Whether any method was used to select for or against, enrich, \
            or screen the material being sequenced. Library Selection method \
            (e.g. random, PCA, cDNA, etc )",
    )
    library_construction: str = Field(
        None,
        description="The name of a library construction approach being used. \
            (eg: '10X v2 sequencing' or 'Smart-seq2')",
    )
    library_preparation: str = Field(
        None,
        description="The general method for sequencing library construction \
            (e.g. KAPA PCR-free).",
    )
    library_level: Optional[str] = Field(
        None, description="Single Cell Sequencing or Bulk Sequencing"
    )
    library_construction_kit_retail_name: str = Field(
        None,
        description="A unique identifier for the kit used to construct a \
            genomic library.  This may include the vendor name, kit name and \
            kit version  (e.g. Agilent sure select Human Exome V8, \
            Twist RefSeq Exome, etc.)",
    )
    library_construction_kit_manufacturer: str = Field(
        None, description="Manufacturer of library construction kit"
    )
    primer: Optional[str] = Field(
        None,
        description="The type of primer used for reverse transcription, e.g. \
        'oligo-dT' or 'random' primer. This allows users to identify content of \
        the cDNA library input e.g. enriched for mRNA.",
    )
    end_bias: Optional[str] = Field(
        None,
        description="The end of the cDNA molecule that is preferentially sequenced, \
            e.g. 3/5 prime tag or end, or the full-length transcript.",
    )
    target_regions: Optional[str] = Field(
        None,
        description="Subset of genes or specific regions of the genome, which are \
            most likely to be involved in the phenotype under study. ",
    )
    rnaseq_strandedness: Optional[str] = Field(
        None,
        description="The strandedness of the library, whether reads come from both \
            strands of the cDNA  or only from the first (antisense) or the second \
            (sense) strand.",
    )
    name: str = Field(
        None,
        description="Name of the library preparation protocol (eg: mRNA-seq \
            library preparation).",
    )
    description: str = Field(
        None,
        description="Description about how a sequencing library was prepared \
            (eg: Library construction method).",
    )
    url: Optional[str] = Field(None, description="A URL to a resource.")
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="One or more attributes that further characterizes this \
            Library Preparation Protocol.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateSequencingProtocol(Protocol):
    """
    Information about the sequencing of a sample.
    """

    sequencing_center: str = Field(
        None, description="Center where sample was sequenced."
    )
    instrument_model: str = Field(
        None,
        description="The name and model of the technology platform used to \
            perform sequencing.",
    )
    read_length: Optional[str] = Field(
        None,
        description="Length of sequencing reads (eg: Long or short or actual number \
            of the read length etc). The number of nucleotides successfully ordered \
                from each side of a nucleic acid fragment  obtained after the \
                    completion of a sequencing process",
    )
    read_pair_number: Optional[str] = Field(
        None,
        description="Denotes whether a submitted FASTQ file contains forward (R1) or \
            reverse (R2) reads for paired-end sequencing. The number that identifies \
            each read direction in a paired-end nucleotide sequencing replications.",
    )
    sequencing_length: Optional[str] = Field(None, description="Long or Short Read.")
    target_coverage: Optional[str] = Field(
        None,
        description="Mean coverage for whole genome sequencing, or mean target \
            coverage for whole exome and targeted sequencing. The number of times a \
            particular locus (site, nucleotide, amplicon, region) was sequenced.",
    )
    reference_annotation: Optional[str] = Field(
        None,
        description="A published genetic sequence that is used as a reference \
            sequence against which other sequences are compared. Reference genome(s) \
            or annotation(s) used for prior analyses (eg: GRCh38.p13).",
    )
    lane_number: Optional[str] = Field(
        None,
        description="he numerical identifier for the lane or machine unit where a \
            sample was located during nucleotide sequencing.",
    )
    flow_cell_id: Optional[str] = Field(
        None,
        description="Flow Cell ID (eg: Experiment ID_Cell 1_Lane_1). The barcode \
            assigned to a flow cell used in nucleotide sequencing.",
    )
    flow_cell_type: Optional[str] = Field(
        None,
        description="Type of flow cell used (e.g. S4, S2 for NovaSeq; PromethION, \
            Flongle for Nanopore). Aparatus in the fluidic subsystem where the \
            sheath and sample meet.  Can be one of several types; jet-in-air, \
            quartz cuvette, or a hybrid of the two.  The sample flows through \
            the center of a fluid column of sheath fluid in the flow cell.",
    )
    umi_barcode_read: Optional[str] = Field(
        None,
        description="The type of read that contains the UMI barcode \
            (Eg: index1/index2/read1/read2).",
    )
    umi_barcode_size: Optional[str] = Field(
        None, description="The size of the UMI identifying barcode (Eg. '10')."
    )
    umi_barcode_offset: Optional[str] = Field(
        None,
        description="The offset in sequence of the UMI identifying barcode. \
            (E.g. '16').",
    )
    cell_barcode_read: Optional[str] = Field(
        None,
        description="The type of read that contains the cell barcode \
            (eg: index1/index2/read1/read2).",
    )
    cell_barcode_offset: Optional[str] = Field(
        None,
        description="The offset in sequence of the cell identifying barcode. \
            (Eg. '0').",
    )
    cell_barcode_size: Optional[str] = Field(
        None, description="The size of the cell identifying barcode (E.g. '16')."
    )
    sample_barcode_read: Optional[str] = Field(
        None,
        description="The type of read that contains the sample barcode \
            (eg: index1/index2/read1/read2).",
    )
    name: Optional[str] = Field(
        None,
        description="Name of the library preparation protocol \
            (eg: mRNA-seq,Whole exome long-read sequencing etc).",
    )
    description: Optional[str] = Field(
        None,
        description="Description about the sequencing protocol \
            (eg: mRNA-seq,Whole exome long-read sequencing etc).",
    )
    url: Optional[str] = Field(None, description="A URL to a resource.")
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="One or more attributes that further characterizes \
            this Sequencing Protocol.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class SequencingProtocol(Protocol):
    """
    Information about the sequencing of a sample.
    """

    sequencing_center: str = Field(
        None, description="Center where sample was sequenced."
    )
    instrument_model: str = Field(
        None,
        description="The name and model of the technology platform used to \
            perform sequencing.",
    )
    read_length: Optional[str] = Field(
        None,
        description="Length of sequencing reads (eg: Long or short or actual number \
            of the read length etc). The number of nucleotides successfully ordered \
                from each side of a nucleic acid fragment  obtained after the \
                    completion of a sequencing process",
    )
    read_pair_number: Optional[str] = Field(
        None,
        description="Denotes whether a submitted FASTQ file contains forward (R1) or \
            reverse (R2) reads for paired-end sequencing. The number that identifies \
            each read direction in a paired-end nucleotide sequencing replications.",
    )
    sequencing_length: Optional[str] = Field(None, description="Long or Short Read.")
    target_coverage: Optional[str] = Field(
        None,
        description="Mean coverage for whole genome sequencing, or mean target \
            coverage for whole exome and targeted sequencing. The number of times a \
            particular locus (site, nucleotide, amplicon, region) was sequenced.",
    )
    reference_annotation: Optional[str] = Field(
        None,
        description="A published genetic sequence that is used as a reference \
            sequence against which other sequences are compared. Reference genome(s) \
            or annotation(s) used for prior analyses (eg: GRCh38.p13).",
    )
    lane_number: Optional[str] = Field(
        None,
        description="he numerical identifier for the lane or machine unit where a \
            sample was located during nucleotide sequencing.",
    )
    flow_cell_id: Optional[str] = Field(
        None,
        description="Flow Cell ID (eg: Experiment ID_Cell 1_Lane_1). The barcode \
            assigned to a flow cell used in nucleotide sequencing.",
    )
    flow_cell_type: Optional[str] = Field(
        None,
        description="Type of flow cell used (e.g. S4, S2 for NovaSeq; PromethION, \
            Flongle for Nanopore). Aparatus in the fluidic subsystem where the \
            sheath and sample meet.  Can be one of several types; jet-in-air, \
            quartz cuvette, or a hybrid of the two.  The sample flows through \
            the center of a fluid column of sheath fluid in the flow cell.",
    )
    umi_barcode_read: Optional[str] = Field(
        None,
        description="The type of read that contains the UMI barcode \
            (Eg: index1/index2/read1/read2).",
    )
    umi_barcode_size: Optional[str] = Field(
        None, description="The size of the UMI identifying barcode (Eg. '10')."
    )
    umi_barcode_offset: Optional[str] = Field(
        None,
        description="The offset in sequence of the UMI identifying barcode. \
            (E.g. '16').",
    )
    cell_barcode_read: Optional[str] = Field(
        None,
        description="The type of read that contains the cell barcode \
            (eg: index1/index2/read1/read2).",
    )
    cell_barcode_offset: Optional[str] = Field(
        None,
        description="The offset in sequence of the cell identifying barcode. \
            (Eg. '0').",
    )
    cell_barcode_size: Optional[str] = Field(
        None, description="The size of the cell identifying barcode (E.g. '16')."
    )
    sample_barcode_read: Optional[str] = Field(
        None,
        description="The type of read that contains the sample barcode \
            (eg: index1/index2/read1/read2).",
    )
    name: Optional[str] = Field(
        None,
        description="Name of the library preparation protocol \
            (eg: mRNA-seq,Whole exome long-read sequencing etc).",
    )
    description: Optional[str] = Field(
        None,
        description="Description about the sequencing protocol \
            (eg: mRNA-seq,Whole exome long-read sequencing etc).",
    )
    url: Optional[str] = Field(None, description="A URL to a resource.")
    has_attribute: Optional[List[Attribute]] = Field(
        None,
        description="One or more attributes that further characterizes \
            this Sequencing Protocol.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateAgent(BaseModel):
    """
    An agent is something that bears some form of responsibility for an
    activity taking place, for the existence of an entity, or for another
    agent's activity. Agents include a Person, Organization, or Software
    that performs an activity.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Agent(BaseModel):
    """
    An agent is something that bears some form of responsibility for an
    activity taking place, for the existence of an entity, or for another
    agent's activity. Agents include a Person, Organization, or Software
    that performs an activity.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateTechnology(BaseModel):
    """
    A Technology is an abstraction that represents the instrument used for an assay.
    The Technology entity captures instrument-specific attributes that are relevant
    for an Experiment entity. The Technology entity may be further characterized by
    its children where each child has fields that are relevant to that particular
    technology.
    """

    alias: Optional[str] = Field(None, description="The alias for an entity.")
    name: Optional[str] = Field(None, description="The name for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Technology(BaseModel):
    """
    A Technology is an abstraction that represents the instrument used for an assay.
    The Technology entity captures instrument-specific attributes that are relevant
    for an Experiment entity. The Technology entity may be further characterized by
    its children where each child has fields that are relevant to that particular
    technology.
    """

    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    name: Optional[str] = Field(None, description="The name for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateWorkflow(BaseModel):
    """
    A Workflow is an abstraction that represents the workflow used to perform an
    analysis. The Workflow entity captures workflow-specific attributes that are
    relevant for an Analysis entity. The Workflow entity may be further characterized
    by its children where each child has fields that are relevant to that
    particular workflow.
    """

    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Workflow(BaseModel):
    """
    A Workflow is an abstraction that represents the workflow used to perform an
    analysis. The Workflow entity captures workflow-specific attributes that are
    relevant for an Analysis entity. The Workflow entity may be further characterized
    by its children where each child has fields that are relevant to that
    particular workflow.
    """

    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateWorkflowStep(BaseModel):
    """
    A Workflow Step represents each individual step performed in a Workflow.
    If the Workflow is a single-step workflow then the Workflow has just one
    Workflow Step entity. If the Workflow is a multi-step workflow then the Workflow
    has a Workflow Step entity for each step. All Workflow step specific attributes
    like parameters, and metadata about execution environment are captured by the
    Workflow Step entity.
    """

    has_parameter: Optional[List[WorkflowParameter]] = Field(
        None,
        description="One or more parameters that are associated with this \
            Workflow Step.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class WorkflowStep(BaseModel):
    """
    A Workflow Step represents each individual step performed in a Workflow.
    If the Workflow is a single-step workflow then the Workflow has just one
    Workflow Step entity. If the Workflow is a multi-step workflow then the Workflow
    has a Workflow Step entity for each step. All Workflow step specific attributes
    like parameters, and metadata about execution environment are captured by the
    Workflow Step entity.
    """

    has_parameter: Optional[List[WorkflowParameter]] = Field(
        None,
        description="One or more parameters that are associated with this \
            Workflow Step.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class WorkflowParameter(BaseModel):
    """
    A key/value pair that represents a parameter used in a Workflow Step.
    """

    key: Optional[str] = Field(
        None, description="Key that represents the parameter name."
    )
    value: Optional[str] = Field(
        None, description="Value corresponding to the the parameter key."
    )


class CreateBiospecimen(BaseModel):
    """
    A Biospecimen is any natural material taken from a biological entity
    (usually a human) for testing, diagnostics, treatment, or research purposes.
    The Biospecimen is linked to the Individual from which the Biospecimen is derived.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    has_individual: Optional[Union[str, CreateIndividual]] = Field(
        None,
        description="The Individual entity from which this Biospecimen was derived.",
    )
    has_anatomical_entity: Optional[Union[str, CreateAnatomicalEntity]] = Field(
        None,
        description="The Anatomical entity, that represents the site, from which the \
            Biospecimen was retrieved. Typically, a concept from \
            Uber-anatomy Ontology (UBERON). For example, 'UBERON:0008307' \
            indicates that the Biospecimen was extracted from the 'Heart Endothelium' \
            of an Individual.",
    )
    has_disease: Optional[Union[List[str], List[CreateDisease]]] = Field(
        None,
        description="The Disease entity that is associated with the Individual. \
            Typically, a concept from Mondo Disease Ontology. For example, \
            'MONDO:0005267' indicates that the Individual suffers from \
            'Heart Disease'.",
    )
    has_phenotypic_feature: Optional[
        Union[List[str], List[CreatePhenotypicFeature]]
    ] = Field(
        None,
        description="The Phenotypic Feature entity that is associated with the \
            Individual. Typically, a concept from Human Phenotype Ontology. \
            For example, 'HP:0100244' indicates that the Individual exhibits \
            'Fibrosarcoma' as one of its phenotype.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Biospecimen(BaseModel):
    """
    A Biospecimen is any natural material taken from a biological entity
    (usually a human) for testing, diagnostics, treatment, or research purposes.
    The Biospecimen is linked to the Individual from which the Biospecimen is derived.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    has_individual: Optional[Union[str, Individual]] = Field(
        None,
        description="The Individual entity from which this Biospecimen was derived.",
    )
    has_anatomical_entity: Optional[Union[str, AnatomicalEntity]] = Field(
        None,
        description="The Anatomical entity, that represents the site, from which the \
            Biospecimen was retrieved. Typically, a concept from \
            Uber-anatomy Ontology (UBERON). For example, 'UBERON:0008307' \
            indicates that the Biospecimen was extracted from the 'Heart Endothelium' \
            of an Individual.",
    )
    has_disease: Optional[Union[List[str], List[Disease]]] = Field(
        None,
        description="The Disease entity that is associated with the Individual. \
            Typically, a concept from Mondo Disease Ontology. For example, \
            'MONDO:0005267' indicates that the Individual suffers from \
            'Heart Disease'.",
    )
    has_phenotypic_feature: Optional[Union[List[str], List[PhenotypicFeature]]] = Field(
        None,
        description="The Phenotypic Feature entity that is associated with the \
            Individual. Typically, a concept from Human Phenotype Ontology. \
            For example, 'HP:0100244' indicates that the Individual exhibits \
            'Fibrosarcoma' as one of its phenotype.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateSample(BaseModel):
    """
    A sample is a limited quantity of something to be used for testing, analysis, \
        inspection, investigation, demonstration, or trial use. A sample is prepared \
        from a Biospecimen (isolate or tissue).
    """

    name: str = Field(
        None,
        description="Name of the sample (eg:GHGAS_Blood_Sample1 or \
            GHGAS_PBMC_RNAseq_S1).",
    )
    description: str = Field(
        None,
        description="Short textual description of the sample \
            (How the sample was collected, sample source,  protocol followed \
            for processing the sample etc).",
    )
    vital_status_at_sampling: Optional[str] = Field(
        None,
        description="Vital Status of an Individual at the point of sampling \
            (eg:'Alive', 'Deceased').",
    )
    tissue: str = Field(
        None,
        description="An anatomical structure consisting of similarly specialized \
            cells and intercellular matrix,  aggregated according to genetically \
            determined spatial relationships, performing a specific function.",
    )
    isolation: Optional[str] = Field(
        None, description="Method or device employed for collecting/isolating a sample."
    )
    storage: Optional[str] = Field(
        None,
        description="Methods by which sample is stored \
            (e.g. frozen in liquid nitrogen).",
    )
    has_individual: Union[str, CreateIndividual] = Field(
        None, description="The Individual from which this Sample was derived from."
    )
    has_anatomical_entity: Optional[Union[str, CreateAnatomicalEntity]] = Field(
        None, description="Anatomical site associated with an entity."
    )

    has_biospecimen: Union[str, Optional[CreateBiospecimen]] = Field(
        None, description="The Biospecimen from which this Sample was prepared from."
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None,
        description="One or more cross-references for this Sample. For example, \
            this Sample may have an EBI BioSamples accession or an \
            EGA Sample accession.",
    )


class Sample(BaseModel):
    """
    A sample is a limited quantity of something to be used for testing, analysis,
    inspection, investigation, demonstration, or trial use. A sample is prepared
    from a Biospecimen (isolate or tissue).
    """

    name: str = Field(
        None,
        description="Name of the sample (eg:GHGAS_Blood_Sample1 or \
            GHGAS_PBMC_RNAseq_S1).",
    )
    description: str = Field(
        None,
        description="Short textual description of the sample \
            (How the sample was collected, sample source,  protocol followed \
            for processing the sample etc).",
    )
    vital_status_at_sampling: Optional[str] = Field(
        None,
        description="Vital Status of an Individual at the point of sampling \
            (eg:'Alive', 'Deceased').",
    )
    tissue: str = Field(
        None,
        description="An anatomical structure consisting of similarly specialized \
            cells and intercellular matrix,  aggregated according to genetically \
            determined spatial relationships, performing a specific function.",
    )
    isolation: Optional[str] = Field(
        None, description="Method or device employed for collecting/isolating a sample."
    )
    storage: Optional[str] = Field(
        None,
        description="Methods by which sample is stored \
            (e.g. frozen in liquid nitrogen).",
    )
    has_individual: Union[str, Individual] = Field(
        None, description="The Individual from which this Sample was derived from."
    )
    has_anatomical_entity: Optional[str] = Field(
        None, description="Anatomical site associated with an entity."
    )

    has_biospecimen: Union[str, Optional[Biospecimen]] = Field(
        None, description="The Biospecimen from which this Sample was prepared from."
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None,
        description="One or more cross-references for this Sample. For example, \
            this Sample may have an EBI BioSamples accession or an \
            EGA Sample accession.",
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateIndividual(BaseModel):
    """
    An Individual is a Person who is participating in a Study.
    """

    gender: Optional[str] = Field(
        None,
        description="Identification as male/masculine, female/feminine or \
            something else,  and association with a (social) role or set of \
            behavioral and cultural traits.",
    )
    sex: BiologicalSexEnum = Field(
        None,
        description="The assemblage of physical properties or qualities by which \
            male is distinguished from female;  the physical difference between \
            male and female; the distinguishing peculiarity of male or female.",
    )
    age: int = Field(None, description="Age of an individual.")
    year_of_birth: Optional[str] = Field(
        None, description="The year in which the individual was born."
    )
    vital_status: VitalStatusEnum = Field(
        None, description="Last known Vital Status of an Individual."
    )
    geographical_region: Optional[str] = Field(
        None,
        description="The geographical region where the Individual is located. \
            Any demarcated area of the Earth; may be determined by both natural \
            and human boundaries.",
    )
    ethnicity: Optional[str] = Field(
        None,
        description="A social group characterized by a distinctive social and \
            cultural tradition  that is maintained from generation to generation.",
    )
    ancestry: Optional[str] = Field(
        None,
        description="A person's descent or lineage, from a person or \
            from a population.",
    )
    has_parent: Optional[Union[List[str], List[CreateIndividual]]] = Field(
        None, description="One or more parent for this Individual."
    )
    has_children: Optional[Union[List[str], List[CreateIndividual]]] = Field(
        None, description="One or more children for this Individual."
    )
    has_disease: Optional[Union[List[str], List[CreateDisease]]] = Field(
        None,
        description="The Disease entity that is associated with this Biospecimen \
            at the time of retrieval from the organism. Typically, a concept from \
            Mondo Disease Ontology. For example, 'MONDO:0003742' indicates that \
            the Individual - from which the Biospecimen was extracted from - suffers \
            from 'Heart Fibrosarcoma'.",
    )
    has_phenotypic_feature: Optional[
        Union[List[str], List[CreatePhenotypicFeature]]
    ] = Field(
        None,
        description="The Phenotypic Feature entity that is associated with this \
            Biospecimen at the time of retrieval from the organism. Typically, a \
            concept from Human Phenotype Ontology. For example, 'HP:0100244' \
            indicates that the Individual - from which the Biospecimen was \
            extracted from - exhibits 'Fibrosarcoma' as one of its phenotype.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    given_name: Optional[str] = Field(None, description="First name.")
    family_name: Optional[str] = Field(None, description="Last name.")
    additional_name: Optional[str] = Field(None, description="Additional name(s).")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Individual(BaseModel):
    """
    An Individual is a Person who is participating in a Study.
    """

    gender: Optional[str] = Field(
        None,
        description="Identification as male/masculine, female/feminine or \
            something else,  and association with a (social) role or set of \
            behavioral and cultural traits.",
    )
    sex: BiologicalSexEnum = Field(
        None,
        description="The assemblage of physical properties or qualities by which \
            male is distinguished from female;  the physical difference between \
            male and female; the distinguishing peculiarity of male or female.",
    )
    age: int = Field(None, description="Age of an individual.")
    year_of_birth: Optional[str] = Field(
        None, description="The year in which the individual was born."
    )
    vital_status: VitalStatusEnum = Field(
        None, description="Last known Vital Status of an Individual."
    )
    geographical_region: Optional[str] = Field(
        None,
        description="The geographical region where the Individual is located. \
            Any demarcated area of the Earth; may be determined by both natural \
            and human boundaries.",
    )
    ethnicity: Optional[str] = Field(
        None,
        description="A social group characterized by a distinctive social and \
            cultural tradition  that is maintained from generation to generation.",
    )
    ancestry: Optional[str] = Field(
        None,
        description="A person's descent or lineage, from a person or \
            from a population.",
    )
    has_parent: Optional[Union[List[str], List[Individual]]] = Field(
        None, description="One or more parent for this Individual."
    )
    has_children: Optional[Union[List[str], List[Individual]]] = Field(
        None, description="One or more children for this Individual."
    )
    has_disease: Optional[Union[List[str], List[Disease]]] = Field(
        None,
        description="The Disease entity that is associated with this Biospecimen \
            at the time of retrieval from the organism. Typically, a concept from \
            Mondo Disease Ontology. For example, 'MONDO:0003742' indicates that \
            the Individual - from which the Biospecimen was extracted from - suffers \
            from 'Heart Fibrosarcoma'.",
    )
    has_phenotypic_feature: Optional[Union[List[str], List[PhenotypicFeature]]] = Field(
        None,
        description="The Phenotypic Feature entity that is associated with this \
            Biospecimen at the time of retrieval from the organism. Typically, a \
            concept from Human Phenotype Ontology. For example, 'HP:0100244' \
            indicates that the Individual - from which the Biospecimen was \
            extracted from - exhibits 'Fibrosarcoma' as one of its phenotype.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    given_name: Optional[str] = Field(None, description="First name.")
    family_name: Optional[str] = Field(None, description="Last name.")
    additional_name: Optional[str] = Field(None, description="Additional name(s).")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class Donor(Individual):
    """
    A Donor is an Individual that participates in a research Study by donating
    a Biospecimen. The use of the Biospecimen is restricted to the consent provided
    by the Donor.
    """

    gender: Optional[str] = Field(
        None,
        description="Identification as male/masculine, female/feminine or something \
            else,  and association with a (social) role or set of behavioral and \
            cultural traits.",
    )
    sex: BiologicalSexEnum = Field(
        None,
        description="The assemblage of physical properties or qualities by which \
            male is distinguished from female;  the physical difference between \
            male and female; the distinguishing peculiarity of male or female.",
    )
    age: int = Field(None, description="Age of an individual.")
    year_of_birth: Optional[str] = Field(
        None, description="The year in which the individual was born."
    )
    vital_status: VitalStatusEnum = Field(
        None,
        description="The state or condition of being living or deceased; also \
            includes the case where the vital status is unknown.",
    )
    geographical_region: Optional[str] = Field(
        None,
        description="The geographical region where the Individual is located. \
            Any demarcated area of the Earth; may be determined by both natural \
            and human boundaries.",
    )
    ethnicity: Optional[str] = Field(
        None,
        description="A social group characterized by a distinctive social and \
            cultural tradition  that is maintained from generation to generation.",
    )
    ancestry: Optional[str] = Field(
        None,
        description="A person's descent or lineage, from a person or from \
            a population.",
    )
    has_parent: Optional[Union[List[str], List[Individual]]] = Field(
        None, description="The parent of an entity."
    )
    has_children: Optional[Union[List[str], List[Individual]]] = Field(
        None, description="The children of an entity."
    )
    has_disease: Optional[Union[List[str], List[Disease]]] = Field(
        None, description="Disease concept that the entity is associated with."
    )
    has_phenotypic_feature: Optional[Union[List[str], List[PhenotypicFeature]]] = Field(
        None,
        description="Phenotypic feature concept that the entity is associated with.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    given_name: Optional[str] = Field(None, description="First name.")
    family_name: Optional[str] = Field(None, description="Last name.")
    additional_name: Optional[str] = Field(None, description="Additional name(s).")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class Population(BaseModel):
    """
    A population is a collection of individuals from the same taxonomic class
    living, counted or sampled at a particular site or in a particular area.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class Family(BaseModel):
    """
    A domestic group, or a number of domestic groups linked through descent
    (demonstrated or stipulated) from a common ancestor, marriage, or adoption.
    """

    has_member: Optional[Union[List[str], List[Individual]]] = Field(
        None,
        description="One or more Individuals that collectively define this Family.",
    )
    has_proband: Optional[Union[str, Individual]] = Field(
        None,
        description="The Individual that is reported to have a disorder which results \
            in the Family being brought into a Study.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    name: Optional[str] = Field(None, description="The name for an entity.")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class Cohort(BaseModel):
    """
    A cohort is a collection of individuals that share a common
    characteristic/observation and have been grouped together for
    a research study/investigation.
    """

    has_member: Optional[Union[List[str], List[Individual]]] = Field(
        None,
        description="One or more Individuals that collectively define this Cohort.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole \
            purpose of referring to that entity in a global scope.",
    )
    name: Optional[str] = Field(None, description="The name for an entity.")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateFile(BaseModel):
    """
    A file represents a data entity that has been generated by a planned process or
    a data transformation activity.
    """

    name: str = Field(None, description="The name for an entity.")
    format: Optional[str] = Field(
        None, description="The format of the file: BAM, SAM, CRAM, BAI, etc."
    )
    size: Optional[str] = Field(None, description="The size of a file in bytes.")
    checksum: Optional[str] = Field(
        None,
        description="A computed value which depends on the contents of a block of \
            data and which is transmitted or  stored along with the data in order \
            to detect corruption of the data.  The receiving system recomputes the \
            checksum based upon the received data and compares this  value with the \
            one sent with the data. If the two values are the same, the receiver has \
            some confidence  that the data was received correctly.",
    )
    file_index: Optional[str] = Field(
        None, description="The index for this file. Commonly for BAM/VCF files."
    )
    category: Optional[str] = Field(
        None,
        description="The category for this file: Whole Genome Sequencing, Whole Exome \
            Sequencing, etc.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class File(BaseModel):
    """
    A file represents a data entity that has been generated by a planned process or
    a data transformation activity.
    """

    name: str = Field(None, description="The name for an entity.")
    format: Optional[str] = Field(
        None, description="The format of the file: BAM, SAM, CRAM, BAI, etc."
    )
    size: Optional[str] = Field(None, description="The size of a file in bytes.")
    checksum: Optional[str] = Field(
        None,
        description="A computed value which depends on the contents of a block of \
            data and which is transmitted or  stored along with the data in order \
            to detect corruption of the data.  The receiving system recomputes the \
            checksum based upon the received data and compares this  value with the \
            one sent with the data. If the two values are the same, the receiver has \
            some confidence  that the data was received correctly.",
    )
    file_index: Optional[str] = Field(
        None, description="The index for this file. Commonly for BAM/VCF files."
    )
    category: Optional[str] = Field(
        None,
        description="The category for this file: Whole Genome Sequencing, Whole Exome \
            Sequencing, etc.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateAnalysis(BaseModel):
    """
    An Analysis is a data transformation that transforms input data to output data.
    The workflow used to achieve this transformation and the individual steps are
    also captured.
    """

    has_input: Optional[Union[List[str], List[CreateFile]]] = Field(
        None, description="The input data File entities used in the Analysis."
    )
    has_study: Optional[Union[str, CreateStudy]] = Field(
        None, description="The Study entity associated with this Analysis."
    )
    has_workflow: Optional[Union[str, CreateWorkflow]] = Field(
        None, description="The Workflow entity associated with this Analysis."
    )
    has_analysis_process: Optional[
        Union[List[str], List[CreateAnalysisProcess]]
    ] = Field(
        None,
        description="One or more Analysis Process entities associated with this \
            Analysis.",
    )
    has_output: Optional[Union[List[str], List[CreateFile]]] = Field(
        None, description="The output data File entities generated by this Analysis."
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    title: Optional[str] = Field(
        None, description="The title that describes an entity."
    )
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Analysis(BaseModel):
    """
    An Analysis is a data transformation that transforms input data to output data.
    The workflow used to achieve this transformation and the individual steps are
    also captured.
    """

    has_input: Optional[Union[List[str], List[File]]] = Field(
        None, description="The input data File entities used in the Analysis."
    )
    has_study: Optional[Union[str, Study]] = Field(
        None, description="The Study entity associated with this Analysis."
    )
    has_workflow: Optional[Union[str, Workflow]] = Field(
        None, description="The Workflow entity associated with this Analysis."
    )
    has_analysis_process: Optional[Union[List[str], List[AnalysisProcess]]] = Field(
        None,
        description="One or more Analysis Process entities associated with this \
            Analysis.",
    )
    has_output: Optional[Union[List[str], List[File]]] = Field(
        None, description="The output data File entities generated by this Analysis."
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    title: Optional[str] = Field(
        None, description="The title that describes an entity."
    )
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateAnalysisProcess(BaseModel):
    """
    An analysis process is a process that describes how one or more Files, from a
    Study, are transformed to another set of Files via a Workflow. The analysis
    process also keeps track of the workflow metadata and the Agent that is running
    the Analysis.
    """

    title: Optional[str] = Field(
        None, description="The title that describes an entity."
    )
    has_input: Optional[Union[List[str], List[CreateFile]]] = Field(
        None, description="The input data File entities used in the Analysis Process."
    )
    has_workflow_step: Optional[CreateWorkflowStep] = Field(
        None,
        description="Workflow Step entity that performs a set of operations on the \
            input data Files to generate a set of output data Files.",
    )
    has_agent: Optional[Union[str, CreateAgent]] = Field(
        None,
        description="The Agent - a software, institution, or human - that is \
            executing or responsible for executing the workflow.",
    )
    has_output: Optional[Union[List[str], List[CreateFile]]] = Field(
        None,
        description="The output data File entities generated by the Analysis Process.",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class AnalysisProcess(BaseModel):
    """
    An analysis process is a process that describes how one or more Files, from a
    Study, are transformed to another set of Files via a Workflow. The analysis
    process also keeps track of the workflow metadata and the Agent that is running
    the Analysis.
    """

    title: Optional[str] = Field(
        None, description="The title that describes an entity."
    )
    has_input: Optional[Union[List[str], List[File]]] = Field(
        None, description="The input data File entities used in the Analysis Process."
    )
    has_workflow_step: Optional[WorkflowStep] = Field(
        None,
        description="Workflow Step entity that performs a set of operations on the \
            input data Files to generate a set of output data Files.",
    )
    has_agent: Optional[Union[str, Agent]] = Field(
        None,
        description="The Agent - a software, institution, or human - that is \
            executing or responsible for executing the workflow.",
    )
    has_output: Optional[Union[List[str], List[File]]] = Field(
        None,
        description="The output data File entities generated by the Analysis Process.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class DataUseCondition(BaseModel):
    """
    Data Use Condition represents the use conditions associated with a Dataset.
    A permission field can have one or more terms that collectively defines the
    data use condition. The modifier determines the interpretation of the
    use permission(s).
    """

    permission: Optional[str] = Field(
        None,
        description="Data use permission. Typically one or more terms from DUO. \
            Preferably descendants of 'DUO:0000001 data use permission'.",
    )
    modifier: Optional[str] = Field(
        None,
        description="Modifier for Data use permission. Preferable descendants of \
            'DUO:0000017 data use modifier'",
    )


class CreateMember(BaseModel):
    """
    Member of an Organization or a Committee.
    """

    email: str = Field(None, description="The email of the Member.")
    telephone: str = Field(None, description="The telephone number of the Member.")
    organization: str = Field(
        None, description="The organization that the Member is part of."
    )
    given_name: Optional[str] = Field(None, description="First name.")
    family_name: Optional[str] = Field(None, description="Last name.")
    additional_name: Optional[str] = Field(None, description="Additional name(s).")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Member(BaseModel):
    """
    Member of an Organization or a Committee.
    """

    email: str = Field(None, description="The email of the Member.")
    telephone: str = Field(None, description="The telephone number of the Member.")
    organization: str = Field(
        None, description="The organization that the Member is part of."
    )
    given_name: Optional[str] = Field(None, description="First name.")
    family_name: Optional[str] = Field(None, description="Last name.")
    additional_name: Optional[str] = Field(None, description="Additional name(s).")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateDataAccessCommittee(BaseModel):
    """
    A group of members that are delegated to grant access to one or more datasets
    after ensuring the minimum criteria for data sharing has been met, and request
    for data use does not raise ethical and/or legal concerns.
    """

    name: str = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    main_contact: Optional[str] = Field(
        None, description="The main contact for the Data Access Committee."
    )
    has_member: Optional[Union[List[str], List[CreateMember]]] = Field(
        None, description="All the members that are part of this Data Access Committee."
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class DataAccessCommittee(BaseModel):
    """
    A group of members that are delegated to grant access to one or more datasets
    after ensuring the minimum criteria for data sharing has been met, and request
    for data use does not raise ethical and/or legal concerns.
    """

    name: str = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    main_contact: Optional[str] = Field(
        None, description="The main contact for the Data Access Committee."
    )
    has_member: Optional[Union[List[str], List[Member]]] = Field(
        None, description="All the members that are part of this Data Access Committee."
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateDataAccessPolicy(BaseModel):
    """
    A Data Access Policy specifies under which circumstances, legal or otherwise,
    a user can have access to one or more Datasets belonging to one or more Studies.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: str = Field(
        None, description="A short description for the Data Access Policy."
    )
    policy_text: str = Field(
        None,
        description="The terms of data use and policy verbiage should be \
            captured here.",
    )
    policy_url: Optional[str] = Field(
        None,
        description="URL for the policy, if available. This is useful if the \
            terms of the policy is made available online at a resolvable URL.",
    )
    has_data_access_committee: Union[str, CreateDataAccessCommittee] = Field(
        None, description="The Data Access Committee linked to this policy."
    )
    has_data_use_condition: Optional[Union[str, List[DataUseCondition]]] = Field(
        None,
        description="Data Use Condition entities that are associated with the \
            Data Access Policy.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class DataAccessPolicy(BaseModel):
    """
    A Data Access Policy specifies under which circumstances, legal or otherwise,
    a user can have access to one or more Datasets belonging to one or more Studies.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: str = Field(
        None, description="A short description for the Data Access Policy."
    )
    policy_text: str = Field(
        None,
        description="The terms of data use and policy verbiage should be \
            captured here.",
    )
    policy_url: Optional[str] = Field(
        None,
        description="URL for the policy, if available. This is useful if the \
            terms of the policy is made available online at a resolvable URL.",
    )
    has_data_access_committee: Union[str, DataAccessCommittee] = Field(
        None, description="The Data Access Committee linked to this policy."
    )
    has_data_use_condition: Optional[Union[str, List[DataUseCondition]]] = Field(
        None,
        description="Data Use Condition entities that are associated with the \
            Data Access Policy.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole \
            purpose of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateDataset(BaseModel):
    """
    A Dataset is a collection of Files that is prepared for distribution.
    """

    has_file: Union[List[str], List[CreateFile]] = Field(
        None,
        description="One or more File entities that collectively are part of \
            this Dataset.",
    )
    has_data_access_policy: Union[str, CreateDataAccessPolicy] = Field(
        None, description="The Data Access Policy that applies to this Dataset."
    )


class DatasetStatusPatch(BaseModel):
    """
    A Dataset is a collection of Files that is prepared for distribution.
    """

    status: Optional[ReleaseStatusEnum] = Field(
        None, description="The status of an entity."
    )


class Dataset(BaseModel):
    """
    A Dataset is a collection of Files that is prepared for distribution.
    """

    title: str = Field(None, description="The title that describes an entity.")
    description: str = Field(None, description="Description of an entity.")
    has_file: Optional[Union[List[str], List[File]]] = Field(
        None,
        description="One or more File entities that collectively are part of \
            this Dataset.",
    )
    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )
    ega_accession: Optional[str] = Field(
        None,
        description="A unique European Genome-Phenome Archive (EGA) identifier \
            assigned to an entity for the sole purpose of referring to that entity \
            within the EGA federated network.",
    )
    has_publication: Optional[Union[List[str], List[Publication]]] = Field(
        None,
        description="One or more Publication entities associated with this Dataset.",
    )
    has_study: Optional[Union[List[str], List[Study]]] = Field(
        None,
        description="One or more Study entities that are referenced by this Dataset.",
    )
    has_analysis: Optional[Union[List[str], List[Analysis]]] = Field(
        None,
        description="One or more Analysis entities that are referenced by \
            this Dataset.",
    )
    has_experiment: Optional[Union[List[str], List[Experiment]]] = Field(
        None,
        description="One or more Experiment entities that are referenced by \
            this Dataset.",
    )
    has_sample: Optional[Union[List[str], List[Sample]]] = Field(
        None,
        description="One or more Sample entities that are referenced by this Dataset.",
    )
    has_data_access_policy: Optional[Union[str, DataAccessPolicy]] = Field(
        None, description="The Data Access Policy that applies to this Dataset."
    )
    status: Optional[ReleaseStatusEnum] = Field(
        None, description="The status of an entity."
    )
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    type: Optional[str] = Field(None, description="The type of Dataset.")
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreatePublication(BaseModel):
    """
    The Publication entity represents a publication. While a publication can be any
    article that is published, the minimum expectation is that the publication
    has a valid DOI.
    """

    title: Optional[str] = Field(None, description="The title for the Publication.")
    abstract: Optional[str] = Field(
        None,
        description="The study abstract that describes the goals.  Can also hold \
            abstract from a publication related to this study",
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="One or more cross-references for this Publication."
    )


class Publication(BaseModel):
    """
    The Publication entity represents a publication. While a publication can be any
    article that is published, the minimum expectation is that the publication
    has a valid DOI.
    """

    title: Optional[str] = Field(None, description="The title for the Publication.")
    abstract: Optional[str] = Field(
        None,
        description="The study abstract that describes the goals.  Can also hold \
            abstract from a publication related to this study",
    )
    id: str = Field(None, description="A PMID or DOI for the Publication.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="One or more cross-references for this Publication."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateAnatomicalEntity(BaseModel):
    """
    Biological entity that is either an individual member of a biological species
    or constitutes the structural organization of an individual member of a
    biological species.
    """

    alias: Optional[str] = Field(None, description="The alias for an entity.")
    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class AnatomicalEntity(BaseModel):
    """
    Biological entity that is either an individual member of a biological species
    or constitutes the structural organization of an individual member of a
    biological species.
    """

    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateDisease(BaseModel):
    """
    A disease is a disposition to undergo pathological processes that exists in an
    organism because of one or more disorders in that organism.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class Disease(BaseModel):
    """
    A disease is a disposition to undergo pathological processes that exists in an
    organism because of one or more disorders in that organism.
    """

    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreatePhenotypicFeature(BaseModel):
    """
    The observable form taken by some character (or group of characters) in an
    individual or an organism, excluding pathology and disease. The detectable
    outward manifestations of a specific genotype.
    """

    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )


class PhenotypicFeature(BaseModel):
    """
    The observable form taken by some character (or group of characters) in an
    individual or an organism, excluding pathology and disease. The detectable
    outward manifestations of a specific genotype.
    """

    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    name: Optional[str] = Field(None, description="The name for an entity.")
    description: Optional[str] = Field(None, description="Description of an entity.")
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class User(BaseModel):
    """
    A user in GHGA.
    """

    email: Optional[str] = Field(None, description="Email of a person.")
    role: Optional[UserRoleEnum] = Field(None, description="The role of the user")
    given_name: Optional[str] = Field(None, description="First name.")
    family_name: Optional[str] = Field(None, description="Last name.")
    additional_name: Optional[str] = Field(None, description="Additional name(s).")
    id: str = Field(
        None, description="An identifier that uniquely represents an entity."
    )
    alias: Optional[str] = Field(None, description="The alias for an entity.")
    xref: Optional[List[str]] = Field(
        None, description="Alternate identifiers for an entity."
    )
    creation_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was created."
    )
    update_date: Optional[str] = Field(
        None, description="Timestamp (in ISO 8601 format) when the entity was updated."
    )


class CreateSubmission(BaseModel):
    """
    A grouping entity that represents information about one or more entities.
    A submission can be considered as a set of inter-related (and inter-connected)
    entities that represent a data submission to GHGA.
    """

    alias: str = Field(None, description="An alias for the submission.")
    has_study: CreateStudy = Field(
        None,
        description="Information about a Study entities associated \
            with this submission.",
    )
    has_project: CreateProject = Field(
        None,
        description="Information about a Project entity associated \
            with this submission.",
    )
    has_sample: List[CreateSample] = Field(
        None,
        description="Information about one or more Sample entities associated \
            with this submission.",
    )
    has_biospecimen: List[CreateBiospecimen] = Field(
        None,
        description="Information about one or more Biospecimen entities associated \
            with this submission.",
    )
    has_individual: List[CreateIndividual] = Field(
        None,
        description="Information about one or more Individual entities associated \
            with this submission.",
    )
    has_experiment: Optional[List[CreateExperiment]] = Field(
        None,
        description="Information about one or more Experiment entities associated \
            with this submission.",
    )
    has_analysis: Optional[List[CreateAnalysis]] = Field(
        None,
        description="Information about one or more Analysis entities associated \
            with this submission.",
    )
    has_file: List[CreateFile] = Field(
        None,
        description="Information about one or more File entities associated \
            with this submission.",
    )

    has_technology: List[Technology] = Field(
        None,
        description="Information about one or more Technology entities associated \
            with this submission.",
    )


class SubmissionStatusPatch(BaseModel):
    """
    An object that can be used to change the status of a Submission.
    """

    status: Optional[SubmissionStatusEnum] = Field(
        None,
        description="The status of a Submission. \
            For example, 'in progress' or 'completed'.",
    )


class Submission(BaseModel):
    """
    A grouping entity that represents information about one or more entities.
    A submission can be considered as a set of inter-related (and inter-connected)
    entities that represent a data submission to GHGA.
    """

    id: str = Field(
        None, description="A internal unique identifier for the Submission."
    )
    has_study: Study = Field(
        None,
        description="Information about a Study entities associated \
            with this submission.",
    )
    has_project: Project = Field(
        None,
        description="Information about a Project entity associated \
            with this submission.",
    )
    has_sample: List[Sample] = Field(
        None,
        description="Information about one or more Sample entities associated \
            with this submission.",
    )
    has_biospecimen: List[Biospecimen] = Field(
        None,
        description="Information about one or more Biospecimen entities associated \
            with this submission.",
    )
    has_individual: List[Individual] = Field(
        None,
        description="Information about one or more Individual entities associated \
            with this submission.",
    )
    has_experiment: Optional[List[Experiment]] = Field(
        None,
        description="Information about one or more Experiment entities associated \
            with this submission.",
    )
    has_analysis: Optional[List[Analysis]] = Field(
        None,
        description="Information about one or more Analysis entities associated \
            with this submission.",
    )
    has_file: List[File] = Field(
        None,
        description="Information about one or more File entities associated \
            with this submission.",
    )
    has_technology: List[Technology] = Field(
        None,
        description="Information about one or more Technology entities associated \
            with this submission.",
    )

    has_data_access_policy: Optional[DataAccessPolicy] = Field(
        None,
        description="The Data Access Policy entity that applies to the data associated \
            with this submission.",
    )
    submission_date: Optional[str] = Field(
        None,
        description="The timestamp (in ISO 8601 format) when submission was \
            marked completed.",
    )
    creation_date: Optional[str] = Field(
        None,
        description="Timestamp (in ISO 8601 format) when the Submission was created.",
    )
    update_date: Optional[str] = Field(
        None,
        description="Timestamp (in ISO 8601 format) when the Submission was updated.",
    )
    status: Optional[SubmissionStatusEnum] = Field(
        None,
        description="The status of a Submission. \
            For example, 'in progress' or 'submitted'.",
    )


class OntologyClassMixin(BaseModel):
    """
    Mixin for entities that represent an class/term/concept from an ontology.
    """

    id: str = Field(
        None,
        description="The Compact UR Identifier (CURIE) that uniquely identifies \
            this ontology class.",
    )
    name: Optional[str] = Field(
        None, description="The name or label (rdfs:label) of an ontology class."
    )
    description: Optional[str] = Field(
        None, description="The description or definition of an ontology class."
    )


class AccessionMixin(BaseModel):
    """
    Mixin for entities that can be assigned a GHGA accession.
    """

    accession: Optional[str] = Field(
        None,
        description="A unique GHGA identifier assigned to an entity for the sole purpose \
            of referring to that entity in a global scope.",
    )


class AttributeMixin(BaseModel):
    """
    Mixin for entities that can have one or more attributes.
    """

    has_attribute: Optional[List[Attribute]] = Field(
        None, description="Key/value pairs corresponding to an entity."
    )


class PublicationMixin(BaseModel):
    """
    Mixin for entities that can have one or more publications.
    """

    has_publication: Optional[str] = Field(
        None, description="The Publication associated with an entity."
    )


class DeprecatedMixin(BaseModel):
    """
    Mixin for entities that can be deprecated.
    """

    replaced_by: Optional[str] = Field(
        None,
        description="Refers to the entity which replaces a currently \
            deprecated entity.",
    )
    deprecation_date: Optional[str] = Field(
        None,
        description="The timestamp (in ISO 8601 format) when the entity \
            was deprecated.",
    )


class StatusMixin(BaseModel):
    """
    Mixin for entities that can have a status.
    """

    status: Optional[ReleaseStatusEnum] = Field(
        None, description="The status of an entity."
    )


class ReleaseMixin(BaseModel):
    """
    Mixin for entities that can be released at a later point in time.
    """

    release_date: Optional[str] = Field(
        None,
        description="The timestamp (in ISO 8601 format) when the entity \
            was released for public consumption.",
    )


Attribute.update_forward_refs()

CreateProject.update_forward_refs()
Project.update_forward_refs()

CreateStudy.update_forward_refs()
Study.update_forward_refs()

CreateExperiment.update_forward_refs()
Experiment.update_forward_refs()

CreateExperimentProcess.update_forward_refs()
ExperimentProcess.update_forward_refs()

CreateLibraryPreparationProtocol.update_forward_refs()
LibraryPreparationProtocol.update_forward_refs()

CreateSequencingProtocol.update_forward_refs()
SequencingProtocol.update_forward_refs()

CreateAgent.update_forward_refs()
Agent.update_forward_refs()

CreateTechnology.update_forward_refs()
Technology.update_forward_refs()

CreateWorkflow.update_forward_refs()
Workflow.update_forward_refs()

CreateWorkflowStep.update_forward_refs()
WorkflowStep.update_forward_refs()

WorkflowParameter.update_forward_refs()

CreateBiospecimen.update_forward_refs()
Biospecimen.update_forward_refs()

CreateSample.update_forward_refs()
Sample.update_forward_refs()

CreateIndividual.update_forward_refs()
Individual.update_forward_refs()

Donor.update_forward_refs()

Population.update_forward_refs()

Family.update_forward_refs()

Cohort.update_forward_refs()

CreateFile.update_forward_refs()
File.update_forward_refs()

CreateAnalysis.update_forward_refs()
Analysis.update_forward_refs()

CreateAnalysisProcess.update_forward_refs()
AnalysisProcess.update_forward_refs()

CreateDataset.update_forward_refs()
Dataset.update_forward_refs()

DataUseCondition.update_forward_refs()

CreateDataAccessPolicy.update_forward_refs()
DataAccessPolicy.update_forward_refs()

CreateDataAccessCommittee.update_forward_refs()
DataAccessCommittee.update_forward_refs()

CreateMember.update_forward_refs()
Member.update_forward_refs()

CreatePublication.update_forward_refs()
Publication.update_forward_refs()

CreateAnatomicalEntity.update_forward_refs()
AnatomicalEntity.update_forward_refs()

CreateDisease.update_forward_refs()
Disease.update_forward_refs()

CreatePhenotypicFeature.update_forward_refs()
PhenotypicFeature.update_forward_refs()

User.update_forward_refs()

CreateSubmission.update_forward_refs()
Submission.update_forward_refs()

OntologyClassMixin.update_forward_refs()

AccessionMixin.update_forward_refs()

AttributeMixin.update_forward_refs()

PublicationMixin.update_forward_refs()

DeprecatedMixin.update_forward_refs()

StatusMixin.update_forward_refs()

ReleaseMixin.update_forward_refs()
