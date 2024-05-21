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

def send_email_bcb(project_name, project_id, email):
    msg = EmailMessage()
    msg["To"] = ['jwang@tome.bio', email] if email else ['jwang@tome.bio']
    msg["From"] = 'bfx@tome.bio'
    msg["Subject"] = f'{project_name} sequencing is finished.'
    msg.set_content(f'Basespace project_id for {project_name} is {project_id}.')

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
    current_run = {}
    while True:
        response = requests.get(f'{bs_api_server}/runs?access_token={bs_access_token}&sortby=DateCreated&SortDir=Desc&limit=5', stream=True)
        for run in response.json().get("Items"):
            samples = {}
            if run["Status"] != "Complete" and run["Status"] != "Failed" and run["Status"] != "Needs Attention":
                current_run[run["V1Pre3Id"]] = run["ExperimentName"]
                print(current_run)
            elif run["V1Pre3Id"] in current_run:
                # store the runinfo and stats
                response = requests.get(
                    f'{bs_api_server}/runs/{run["V1Pre3Id"]}/sequencingstats?access_token={bs_access_token}', stream=True)
                run_json = {"bsrunid": run["V1Pre3Id"], "q30_percentage": format(response.json().get("PercentGtQ30"), ".2f")}

                response = requests.get(
                    f'{bs_api_server}/datasets?InputRuns={run["V1Pre3Id"]}&access_token={bs_access_token}&limit=2048', stream=True)
                for item in response.json().get("Items"):
                    project = item.get("Project").get("Name")
                    if project != "Unindexed Reads":
                        samples[project] = item.get("Project").get("Id")

                send_email(run["ExperimentName"], samples.keys())
                subprocess.call("python /home/ubuntu/bin/tbOnT/update_SRTB.py --srtb %s" % (run["V1Pre3Id"]), shell=True)
                del current_run[run["V1Pre3Id"]]

                benchling = Benchling(url=api_url, auth_method=ApiKeyAuth(api_key))
                for s, id in samples.items():

                    # Query BTB or CTB entity
                    ngs_id = re.sub(".*([C|B]TB\d+).*", "\\1", s)
                    entity = benchling.custom_entities.list(name=ngs_id)
                    ngs_name = entity.first()

                    # change status of entity to sequencing complete or data analysis
                    if "BTB" in s:
                        update = CustomEntityUpdate(fields=fields({"job status": {"value": "sfso_6aKzgWvN"}}))
                        benchling.custom_entities.update(entity_id=ngs_name.id, entity=update)
                    elif "CTB" in s:
                        update = CustomEntityUpdate(fields=fields({"Status": {"value": "sfso_3jyKGSJ6"}}))
                        benchling.custom_entities.update(entity_id=ngs_name.id, entity=update)

                        # send email to the BCB person
                        bcb_email = ngs_name.fields.get("Bioinformatician (email)").text_value
                        send_email_bcb(s, id, bcb_email)

                    # only process BTB project
                    if "BTB" not in s:
                        continue

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
                    if ngs_name.fields.get("NGS assay").text_value == "rhAmpSeq":
                        subprocess.call("python /home/ubuntu/bin/tbOnT/tbrhAmpSeq.py -m %s -i %s -p 8 -o %s" % (
                            s, s, pipeline_run_name), shell=True)
                    else:
                        subprocess.call("python /home/ubuntu/bin/tbOnT/tbAmpSeq.py -m %s -i %s -p 8 -o %s" % (
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
                                                               "AmpSeq Project Name": {"value": ngs_name.id},
                                                               "run start": {"value": run_json["run start"]},
                                                               "run end": {"value": run_json["run end"]},
                                                               "run status": {"value": "Complete"}}))
                    pipeline_run_entity = benchling.custom_entities.create(entity).id

                    ## check the accessibility of the ELN
                    try:
                        requests = benchling.requests.get_by_id(run_json["project_name"])
                    except:
                        pass
                    else:
                        update = CustomEntityUpdate(fields=fields({"ELN entry": {"value": run_json["project_name"]}}))
                        benchling.custom_entities.update(entity_id=pipeline_run_entity, entity=update)

                    # push the data to quilt
                    if ngs_name.fields.get("NGS assay").text_value == "rhAmpSeq":
                        p = subprocess.Popen(
                            "python /home/ubuntu/bin/tbOnT/tbrhAmpSeq.quilt.py -m %s -i %s" % (pipeline_run_name, pipeline_run_name),
                            stdout=subprocess.PIPE, shell=True)
                    else:
                        p = subprocess.Popen(
                            "python /home/ubuntu/bin/tbOnT/quilt.py -m %s -i %s" % (pipeline_run_name, pipeline_run_name),
                            stdout=subprocess.PIPE, shell=True)
                    quilt_link = p.communicate()[0].decode('utf-8').rstrip()
                    if quilt_link:
                        update = CustomEntityUpdate(fields=fields({"analysis result URL link": {"value": quilt_link},
                                                                   "job status": {"value": "sfso_NSDzT3ki"}}))
                        benchling.custom_entities.update(entity_id=ngs_name.id, entity=update)

                    # backup the data to S3
                    subprocess.call("aws s3 sync %s s3://tb-ngs-raw/MiSeq/%s --quiet " % (s, s), shell=True)
                    subprocess.call("aws s3 sync %s s3://tb-ngs-analysis/%s --quiet" % (pipeline_run_name, pipeline_run_name), shell=True)

        time.sleep(3600)
