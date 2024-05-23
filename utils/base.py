import psycopg2
import os
from dotenv import load_dotenv

load_dotenv()

username = os.getenv('WAREHOUSE_USERNAME')
password = os.getenv('WAREHOUSE_PASSWORD')
url = os.getenv('WAREHOUSE_URL')

conn = psycopg2.connect(dbname="warehouse", user=username, password=password, port=5432, sslmode='verify-ca', host=url, connect_timeout=0)
cur = conn.cursor()

api_key = os.getenv('API_KEY')
api_url = os.getenv('API_URL')

registry_id = os.getenv('REGISTRY_ID')

schema_id = os.getenv('AMPSEQ_PIPELINE_RUN_SCHEMA_ID')
folder_id = os.getenv('AMPSEQ_PIPELINE_RUN_FOLDER_ID')

srtb_schema_id = os.getenv('RUN_SCHEMA_ID')
srtb_folder_id = os.getenv('RUN_FOLDER_ID')

result_schema_id = os.getenv('AMPSEQ_RESULTS_SCHEMA_ID')
result_project_id = os.getenv('AMPSEQ_RESULTS_PROJECT_ID')

bs_api_server = os.getenv('apiServer')
bs_access_token = os.getenv('accessToken')

aws_ses_id = os.getenv('AWS_SES_ID')
aws_ses_password = os.getenv('AWS_SES_PASSWORD')