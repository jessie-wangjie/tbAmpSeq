# basespace tracking

import requests
import time
import smtplib
from utils.base import *


def send_email(run_id, samples):
    sender = 'jwang@tome.bio'
    receivers = ['jwang@tome.bio']

    subject = run_id + " is finished."
    message = "Subject:" + subject + "\n" + "\n".join(samples)

    try:
        server = smtplib.SMTP('email-smtp.us-east-1.amazonaws.com', 587)
        server.ehlo()
        server.starttls()
        server.login(aws_ses_id, aws_ses_password)
        server.sendmail(sender, receivers, message)
        print("Successfully sent email")
    except Exception as exception:
        print("Error: %s!\n\n" % exception)


def main():
    current_run = {"TB_MISEQ_000071": "https://api.basespace.illumina.com/v2/runs/247010831"}
    while True:
        response = requests.get(
            f'{bs_api_server}/runs?access_token={bs_access_token}&sortby=DateCreated&SortDir=Desc&limit=5', stream=True)
        for run in response.json().get("Items"):
            samples = {}
            if run["Status"] != "Complete":
                current_run[run["ExperimentName"]] = run["Href"]
            elif run["ExperimentName"] in current_run:
                response = requests.get(
                    f'{current_run[run["ExperimentName"]]}/properties/Input.BioSamples/items?access_token={bs_access_token}&limit=1000',
                    stream=True)
                for item in response.json().get("Items"):
                    project = item.get("BioSample").get("DefaultProject").get("Name")
                    samples[project] = item.get("BioSample").get("DefaultProject").get("Id")
                # send_email(run["ExperimentName"], samples.keys())
                del current_run[run["ExperimentName"]]
                return samples

        time.sleep(18000)
