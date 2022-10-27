import psycopg2
import os
from dotenv import load_dotenv, dotenv_values

config = dotenv_values(".test.env")
print(config.keys())

username = config['WAREHOUSE_USERNAME']
password = config['WAREHOUSE_PASSWORD']
url = config['WAREHOUSE_URL']

conn = psycopg2.connect(f"dbname=warehouse user={username} password={password} port=5432 host={url}")
cur = conn.cursor()

api_key = config['API_KEY']
api_url = config['API_URL']

schema_id = config['AMPSEQ_PIPELINE_RUN_SCHEMA_ID']
folder_id = config['AMPSEQ_PIPELINE_RUN_FOLDER_ID']
registry_id = config['AMPSEQ_PIPELINE_RUN_REGISTRY_ID']

result_schema_id = config['AMPSEQ_RESULTS_SCHEMA_ID']
result_project_id = config['AMPSEQ_RESULTS_PROJECT_ID']
