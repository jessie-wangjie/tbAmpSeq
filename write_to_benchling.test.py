#!/usr/bin/env python
"""
On-Target analysis for Amp-seq and 3Primer-Seq
Information from Benchling
"""

import argparse
import glob
import json
import os

import pandas as pd
from benchling_api_client.models.naming_strategy import NamingStrategy
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.benchling import Benchling
from benchling_sdk.helpers.serialization_helpers import fields
from benchling_sdk.models import CustomEntityCreate, AssayResultCreate, AssayFieldsCreate

from utils.base_test import *


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Write JSON to Benchling',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')

    args = parser.parse_args()
    tbid = args.m

    # read NGS information
    if os.path.exists(os.path.join(tbid, tbid + ".run.json")):
        ngs_stats = pd.read_json(os.path.join(tbid, tbid + ".run.json"), typ="series")

    # create pipeline run entity
    # TODO: to check run suffix
    benchling = Benchling(url=api_url, auth_method=ApiKeyAuth(api_key))

    # check if the entity exists
    entity = benchling.custom_entities.list(name=tbid + "a")
    if entity.estimated_count > 0:
        pipeline_run_entity = entity.first().id
    else:
        entity = CustomEntityCreate(schema_id=schema_id, folder_id=folder_id, registry_id=registry_id,
                                    naming_strategy=NamingStrategy.NEW_IDS, name=tbid + "a",
                                    fields=fields(
                                        {"Genomics AmpSeq Project Queue": {"value": tbid},
                                         "pipeline Name": {"value": "tbAmpseq"},
                                         "github address": {"value": "https://github.com/tomebio/tbOnT"},
                                         "ELN entry": {"value": "etr_D1dGWLNY"},
                                         "AmpSeq Project Name": {"value": "bfi_YfJSzglk"},
                                         "run start": {"value": ngs_stats["run start"]},
                                         "run end": {"value": ngs_stats["run end"]},
                                         "run status": {"value": "Complete"}}))
        pipeline_run_entity = benchling.custom_entities.create(entity)

    # insert the cs2 stats to benchling
    files = glob.glob(tbid + "/*/CRISPResso_stats.json")
    for s in files:
        data = json.load(open(s))
        del data["well"]
        del data["plate"]
        del data["email"]
        data["ampseq_pipeline_run"] = pipeline_run_entity.id
        data["aaan_id"] = "bfi_bBH1gNMA"
        data["ngs_tracking"] = "bfi_YfJSzglk"
        data["ppid"] = "bfi_hoqIh4uj"
        data["project_name"] = ["etr_D1dGWLNY"]
        row = AssayResultCreate(schema_id=result_schema_id, fields=AssayFieldsCreate.from_dict(data), project_id=result_project_id)
        benchling.assay_results.create([row])

if __name__ == "__main__":
    main()
