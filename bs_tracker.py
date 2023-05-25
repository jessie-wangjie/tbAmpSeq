# basespace tracking
import re
import smtplib
import subprocess
import time
from email.message import EmailMessage
import pandas as pd
import requests
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.benchling import Benchling
from benchling_sdk.helpers.serialization_helpers import fields
from benchling_sdk.models import CustomEntityUpdate, CustomEntityCreate
from benchling_api_client.models.naming_strategy import NamingStrategy
from utils.base import *


def send_email(run_id, samples):
    msg = EmailMessage()
    msg["From"] = 'bfx@tome.bio'
    msg["To"] = ['jwang@tome.bio']
    msg["Subject"] = f'{run_id} is finished.'
    msg.set_content("Projects in " + run_id + ":\n\n" + "\n".join(samples))

    try:
        server = smtplib.SMTP('email-smtp.us-east-1.amazonaws.com', 587)
        server.ehlo()
        server.starttls()
        server.login(aws_ses_id, aws_ses_password)
        server.send_message(msg)
        print("Successfully sent email")
    except Exception as exception:
        print("Error: %s!\n\n" % exception)


if __name__ == '__main__':
    current_run = {"TB_MISEQ_000159":259509263}
    while True:
        response = requests.get(
            f'{bs_api_server}/runs?access_token={bs_access_token}&sortby=DateCreated&SortDir=Desc&limit=20', stream=True)
        for run in response.json().get("Items"):
            samples = {}
            if run["Status"] != "Complete":
                current_run[run["ExperimentName"]] = run["V1Pre3Id"]
                print(current_run)
            elif run["ExperimentName"] in current_run:
                print(current_run[run["ExperimentName"]])
                # store the runinfo and stats
                response = requests.get(
                    f'{bs_api_server}/runs/{current_run[run["ExperimentName"]]}/sequencingstats?access_token={bs_access_token}',
                    stream=True)
                run_json = {"bsrunid": run["ExperimentName"], "q30_percentage": format(response.json().get("PercentGtQ30"), ".2f")}

                response = requests.get(
                    f'{bs_api_server}/datasets?InputRuns={current_run[run["ExperimentName"]]}&access_token={bs_access_token}&limit=1000',
                    stream=True)
                for item in response.json().get("Items"):
                    project = item.get("Project").get("Name")
                    if project != "Unindexed Reads":
                        samples[project] = item.get("Project").get("Id")

                send_email(run["ExperimentName"], samples.keys())
                del current_run[run["ExperimentName"]]

                print(samples.items())
                benchling = Benchling(url=api_url, auth_method=ApiKeyAuth(api_key))
                for s, id in samples.items():
                    # check if it's Ampseq data
                    if "BTB" not in s:
                        continue

                    # change status of NGS tracking entity to sequencing complete
                    ngs_id = re.sub(".*(BTB\d+).*", "\\1", s)
                    entity = benchling.custom_entities.list(name=ngs_id)
                    ngs_name = entity.first()
                    update = CustomEntityUpdate(fields=fields({"job status": {"value": "sfso_6aKzgWvN"}}))
                    updated_entity = benchling.custom_entities.update(entity_id=ngs_name.id, entity=update)

                    # check if the pipeline result entity exists
                    entity = benchling.custom_entities.list(name_includes=s, sort="name:desc")
                    if entity.estimated_count > 0:
                        pipeline_run_entity = entity.first()
                        pipeline_run_name = pipeline_run_entity.name[:-1] + chr(ord(pipeline_run_entity.name[-1])+1)
                    else:
                        pipeline_run_name = s + "a"

                    # download fastq files from basespace
                    subprocess.call("bs download project -i %s -o %s --extension=fastq.gz" % (id, s), shell=True)

                    # run CRISPresso2
                    os.makedirs(pipeline_run_name, exist_ok=True)
                    pd.Series(run_json).to_json(os.path.join(pipeline_run_name, "run.json"))
                    subprocess.call("python /home/ubuntu/bin/tbOnT/tbAmpSeq.B2B.coordinates.py -m %s -i %s -p 8 -o %s" % (
                    s, s, pipeline_run_name), shell=True)
                    run_json = pd.read_json(os.path.join(pipeline_run_name, "run.json"), typ="series")

                    # get tbAmpseq commit id
                    p = subprocess.Popen("git -C /data/bin/tbOnT/ rev-parse HEAD", stdout=subprocess.PIPE, shell=True)
                    commit = p.communicate()[0].decode('utf-8').rstrip()

                    entity = CustomEntityCreate(schema_id=schema_id, folder_id=folder_id, registry_id=registry_id,
                                                naming_strategy=NamingStrategy.NEW_IDS, name=pipeline_run_name,
                                                fields=fields({"Genomics AmpSeq Project Queue": {"value": s},
                                                               "pipeline Name": {"value": "tbAmpseq"},
                                                               "github address": {"value": "https://github.com/tomebio/tbOnT"},
                                                               "git commit": {"value": commit},
                                                               "ELN entry": {"value": run_json["project_name"]},
                                                               "AmpSeq Project Name": {"value": ngs_name.id},
                                                               "run start": {"value": run_json["run start"]},
                                                               "run end": {"value": run_json["run end"]},
                                                               "run status": {"value": "Complete"}}))
                    pipeline_run_entity = benchling.custom_entities.create(entity).id

                    # push the data to quilt
                    p = subprocess.Popen(
                        "python /home/ubuntu/bin/tbOnT/quilt.py -m %s -i %s" % (pipeline_run_name, pipeline_run_name),
                        stdout=subprocess.PIPE, shell=True)
                    quilt_link = p.communicate()[0].decode('utf-8').rstrip()
                    if quilt_link:
                        update = CustomEntityUpdate(fields=fields({"analysis result URL link": {"value": quilt_link},
                                                                   "job status": {"value": "sfso_NSDzT3ki"}}))
                        updated_entity = benchling.custom_entities.update(entity_id=ngs_name.id, entity=update)

                    # backup the data to S3
                    subprocess.call("aws s3 --profile=jwang sync %s s3://tb-ngs-raw/MiSeq/%s --quiet " % (s, s), shell=True)
                    # subprocess.call("aws s3 --profile=jwang sync %s s3://tb-ngs-quilt/%s/fastq/ --quiet" % (s, ngs_id), shell=True)
                    subprocess.call("aws s3 --profile=jwang sync %s s3://tb-ngs-analysis/%s --quiet" % (pipeline_run_name, s), shell=True)

        time.sleep(3600)
