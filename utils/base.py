import psycopg2
import os
from dotenv import load_dotenv

load_dotenv()

username = os.getenv('WAREHOUSE_USERNAME')
password = os.getenv('WAREHOUSE_PASSWORD')
url = os.getenv('WAREHOUSE_URL')

# conn = psycopg2.connect(f"dbname=warehouse user={username} password={password} port=5432 host={url}")
# cur = conn.cursor()

api_key = os.getenv('API_KEY')
api_url = os.getenv('API_URL')

schema_id = os.getenv('AMPSEQ_PIPELINE_RUN_SCHEMA_ID')
folder_id = os.getenv('AMPSEQ_PIPELINE_RUN_FOLDER_ID')
registry_id = os.getenv('AMPSEQ_PIPELINE_RUN_REGISTRY_ID')

# test
test_username = os.getenv('TEST_WAREHOUSE_USERNAME')
test_password = os.getenv('TEST_WAREHOUSE_PASSWORD')
test_url = os.getenv('TEST_WAREHOUSE_URL')

conn = psycopg2.connect(f"dbname=warehouse user={test_username} password={test_password} port=5432 host={test_url}")
cur = conn.cursor()

test_api_key = os.getenv('TEST_API_KEY')
test_api_url = os.getenv('TEST_API_URL')
