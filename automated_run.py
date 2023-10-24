import logging
import requests

from bs_tracker import send_email, invoke_pipeline
from utils.base import *

logger = logging.getLogger()
logger.setLevel(logging.INFO)

if __name__ == '__main__':
    run_id = os.environ['RUN_ID']
    project_id = os.environ['PROJECT_ID']
    project_name = os.environ['PROJECT_NAME']
    q30 = os.environ['Q30']
    send_email(os.environ['EXPERIMENT_NAME'], project_name)
    logger.info(f'Running {run_id}')

    invoke_pipeline(sample_list={project_name: project_id}, run_info_json=run_json)
