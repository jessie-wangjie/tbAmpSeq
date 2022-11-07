# basespace tracking

import requests
import time
import smtplib
from utils.base import *


def send_email(run_id, samples):
    sender = 'jessie.wangjie@gmail.com'
    receivers = ['jwang@tome.bio']

    subject = run_id + " is finished."
    message = "Subject:" + subject + "\n" + "\n".join(samples)

    try:
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.ehlo()
        server.starttls()
        server.login("jessie.wangjie@gmail.com", "wquyfnwwpbygpxmt")
        server.sendmail(sender, receivers, message)
        print("Successfully sent email")
    except Exception as exception:
        print("Error: %s!\n\n" % exception)


if __name__ == '__main__':
    current_run = {"TB_MISEQ_000071":"https://api.basespace.illumina.com/v2/runs/247010831"}
    while True:
        response = requests.get(
            f'{bs_api_server}/runs?access_token={bs_access_token}&sortby=DateCreated&SortDir=Desc&limit=5', stream=True)
        for run in response.json().get("Items"):
            samples = {}
            if run["Status"] != "Complete":
                current_run[run["ExperimentName"]] = run["Href"]
            elif run["ExperimentName"] in current_run:
                response = requests.get(
                    f'{current_run[run["ExperimentName"]]}/properties/Input.BioSamples/items&limit=1000?access_token={bs_access_token}',
                    stream=True)
                for item in response.json().get("Items"):
                    project = item.get("BioSample").get("DefaultProject").get("Name")
                    if project not in samples:
                        samples[project] = [item.get("BioSample").get("BioSampleName")]
                    else:
                        samples[project].append(item.get("BioSample").get("BioSampleName"))
                    print(project)
                send_email(run["ExperimentName"], samples.keys())
                del current_run[run["ExperimentName"]]

        print(current_run)
        time.sleep(18000)
