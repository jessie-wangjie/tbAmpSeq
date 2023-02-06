# basespace tracking
import subprocess
import requests
import time
import smtplib
from email.message import EmailMessage
import pandas as pd
from utils.base import *


def send_email(run_id, samples):
    msg = EmailMessage()
    msg["From"] = 'jwang@tome.bio'
    msg["To"] = ['jwang@tome.bio', 'wwang@tome.bio']
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
    current_run = {"TB_MISEQ_000104": "252285035"}
    while True:
        response = requests.get(
            f'{bs_api_server}/runs?access_token={bs_access_token}&sortby=DateCreated&SortDir=Desc&limit=5', stream=True)
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

                # download fastq files from basespace
                print(samples.items())
                for s, id in samples.items():
                    subprocess.call("bs download project -i %s -o %s --extension=fastq.gz" % (id, s), shell=True)
                    os.makedirs(os.path.join(s + "_tbAmpSeq"), exist_ok=True)
                    pd.Series(run_json).to_json(os.path.join(s + "_tbAmpSeq", s + ".run.json"))
                    subprocess.call("python /home/ubuntu/bin/tbOnT/tbAmpSeq.B2B.coordinates.py -m %s -i %s -p 8 -o %s" % (s, s, s + "_tbAmpSeq"), shell=True)

        time.sleep(7200)
