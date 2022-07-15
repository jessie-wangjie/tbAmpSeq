import psycopg2
import os
from dotenv import load_dotenv

load_dotenv()

username = os.getenv('WAREHOUSE_USERNAME')
password = os.getenv('WAREHOUSE_PASSWORD')
url = os.getenv('WAREHOUSE_URL')

conn = psycopg2.connect(f"dbname=warehouse user={username} password={password} port=5432 host={url}")
cur = conn.cursor()
