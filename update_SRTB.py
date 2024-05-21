import argparse
import requests
import re
from utils.base import *
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.models import CustomEntityUpdate, CustomEntityCreate
from benchling_sdk.helpers.serialization_helpers import fields
from benchling_api_client.models.naming_strategy import NamingStrategy


def update_srtb(id):

    basic_keys = ['Name', 'ExperimentName', 'Status', 'DateCreated', 'DateModified', 'FlowcellBarcode', 'ReagentBarcode']
    instrument_keys = ['Name', 'Type']
    stats_keys = ['Chemistry', 'ErrorRate', 'ErrorRateR1', 'ErrorRateR2', 'IntensityCycle1', 'IsIndexed', 'MaxCycleCalled',
                  'MaxCycleExtracted', 'MaxCycleScored', 'MinCycleCalled', 'MinCycleExtracted', 'MinCycleScored',
                  'NonIndexedErrorRate', 'NonIndexedIntensityCycle1', 'NonIndexedPercentAligned', 'NonIndexedPercentGtQ30',
                  'NonIndexedProjectedTotalYield', 'NonIndexedYieldTotal', 'NumCyclesIndex1', 'NumCyclesIndex2',
                  'NumCyclesRead1', 'NumCyclesRead2', 'NumLanes', 'NumReads', 'NumSurfaces', 'NumSwathsPerLane',
                  'NumTilesPerSwath', 'PercentAligned', 'PercentGtQ30', 'PercentGtQ30R1', 'PercentGtQ30R2',
                  'PercentGtQ30Last10Cycles', 'PercentPf', 'PercentResynthesis', 'PhasingR1', 'PhasingR2', 'PrePhasingR1',
                  'PrePhasingR2', 'ProjectedTotalYield', 'ReadsPfTotal', 'ReadsTotal', 'YieldTotal', 'Clusters',
                  'ClustersPf', 'ClusterDensity', 'Occupancy']

    benchling_keys = ['NumCyclesRead1', 'NumCyclesIndex1', 'NumCyclesIndex2', 'NumCyclesRead2', 'PercentGtQ30', 'PercentPf',
                      'PercentAligned', 'YieldTotal', 'ReadsTotal', 'ReadsPfTotal', 'NumLanes', 'Clusters', 'ClustersPf',
                      'ClusterDensity', 'PercentLoadingConcentration']

    # run link
    json_entity = {"Link": {"value": "https://basespace.illumina.com/run/" + id},
                   "V1Pre3Id": {"value": int(id)}, "Projects": {"value": []}}

    # get run basic info
    response = requests.get(f'{bs_api_server}/runs/{id}?access_token={bs_access_token}', stream=True)
    run = response.json()

    for k in basic_keys:
        json_entity[k] = {"value": run.get(k, '')}

    # get run instrument info
    for k in instrument_keys:
        json_entity["Instrument" + k] = {"value": run["Instrument"].get(k, '')}

    # get run stats
    stats = requests.get(f'{bs_api_server}/runs/{id}/sequencingstats?access_token={bs_access_token}', stream=True)
    for k in benchling_keys:
        json_entity[k] = {"value": str(stats.json().get(k, ''))}

    # get projects in each run
    dataset = requests.get(f'{bs_api_server}/datasets?InputRuns={id}&access_token={bs_access_token}&limit=2048', stream=True)
    samples = {}
    for item in dataset.json().get("Items"):
        project = item.get("Project").get("Name")
        if project != "Unindexed Reads":
            samples[project] = item.get("Project").get("Id")
    json_entity["ProjectNames"] = {"value": ",".join([k + ":" + v for k, v in samples.items()])}

    # get benchling id for each project
    run_project = {}
    benchling = Benchling(url=api_url, auth_method=ApiKeyAuth(api_key))
    for s, id in samples.items():
        # change status of NGS tracking entity to sequencing complete
        tb_name = re.sub(".*([C|B]TB\d+).*", "\\1", s)
        entity = benchling.custom_entities.list(name=tb_name)
        if entity.estimated_count > 0:
            tb_entity = entity.first()
            json_entity["Projects"]["value"].append(tb_entity.id)
            run_project[tb_entity.id] = id

    # create or update run entity
    entity_name = run["ExperimentName"] + "_" + run["V1Pre3Id"]
    run_entities = benchling.custom_entities.list(name=entity_name)
    if run_entities.estimated_count > 0:
        run_entity = run_entities.first()
        update_fields = CustomEntityUpdate(fields=fields(json_entity))
        benchling.custom_entities.update(entity_id=run_entity.id, entity=update_fields)
    else:
        entity = CustomEntityCreate(schema_id=schema_id, folder_id=folder_id, registry_id=registry_id,
                                    naming_strategy=NamingStrategy.NEW_IDS, name=entity_name,
                                    fields=fields(json_entity))
        run_entity = benchling.custom_entities.create(entity)

    # update the run info for the project entities, only for CTB
    for id in json_entity["Projects"]["value"]:
        project_entity = benchling.custom_entities.get_by_id(id)
        if "CTB" not in project_entity.entity_registry_id:
            continue

        # add run
        run_set = set(project_entity.fields.get("Sequencing RUN ID").value + [run_entity.id])
        # add project id
        if project_entity.fields.get("Sequencing Project ID").value:
            project_set = set(project_entity.fields.get("Sequencing Project ID").value.split(",") + [run_project[id]])
        else:
            project_set = {run_project[id]}

        update_fields = CustomEntityUpdate(fields=fields({"Sequencing RUN ID": {"value": list(run_set)},
                                                          "Sequencing Project ID": {"value": ",".join(project_set)}}))
        benchling.custom_entities.update(entity_id=id, entity=update_fields)


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Create or update SRTB")

    # Add the arguments
    parser.add_argument("--srtb", required=True, type=str)

    # Parse the arguments
    args = parser.parse_args()

    update_srtb(args.srtb)