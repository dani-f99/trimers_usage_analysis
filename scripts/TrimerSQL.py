
# Imports  
# if needed python3 -m pip install PyMySQL[rsa]
import pymysql
from os.path import exists
import pandas as pd
import os


# Function that conneccting to the mySQL databse and downloading the required
class mysql():
    # Initiation of the SQL server connection
    def __init__(self,
                input_host: str,
                input_port: int,
                input_username: str,
                input_password: str,
                input_database: str):
        
        """
        input_host : str -> sql server host ip, in string format.
        input_port : str -> sql server connection port, in string format.
        input_username : str -> sql server conenction username for connection.
        input_password : str -> sql server connection password for connection.
        input_database : str -> sql server database (from which to download the tables)
        """

        self.db_info = [input_host, input_port, input_username, input_password, input_database]
        
        # Connecting to the mySQL server
        try:
            self.conn = pymysql.connect(host = self.db_info[0],
                                        port = self.db_info[1],
                                        user = self.db_info[2],
                                        password = self.db_info[3],
                                        db = self.db_info[4])
            
            if isinstance(self.conn, pymysql.connections.Connection):
                conn_status = "Successfully connected to MySQL server."
                print(conn_status)
            else:
                conn_status = "Error while connecnting to mysql server, cheek config.ini..."
                raise TypeError(conn_status)
            
        except:
            conn_status = "Error while connecnting to mysql server, cheek function inputs..."
            raise TypeError(conn_status)
    
    # Creating folder and importing the database
    def import_tables(self,
                      tables = ["sample_metadata", "sequences", "sequence_collapse"]):
        
        """
        tabels : list of strings -> names of the required tables from the sql server, will be downloaded in csv format.
        """

        # Cheeking if the imports folder exists
        self.raw_path = "trimers_raw_tables"
        raw_exists = exists(self.raw_path)

        if exists(self.raw_path):
            print("raw dataset folder already exists, continuing...")
        else:
            print("raw dataset folder doesn't exists, creating folder...")
            os.makedirs(self.raw_path)

            if exists(self.raw_path) == True:
                print("raw dataset folder created successfully, continuing...")
            else:
                raise TypeError("error while creating folder, halting...")

        # Importing the tables from the mySQL server    
        self.tables = tables
        len_tables = len(self.tables)
        path_csv = self.raw_path

        num = 1
        for i in self.tables:

            if i in ["sample_metadata", "sequence_collapse"]:
                qry = "SELECT * FROM {t_name};".format(t_name=i)
            else:
                qry = "SELECT * FROM {t_name} WHERE functional = 1 AND sample_id IS NOT NULL AND clone_id IS NOT NULL;".format(t_name=i)
            

            print("This process  may take a while. \n")
            table_name = "{db}.{table}".format(db=self.db_info[4], table=i)
            temp_count = "\n ({n}/{len})".format(n=num, len=len_tables)
            print(table_name)
            
            if exists(path_csv + "\\{t_name}.csv".format(t_name=table_name)):
                print(table_name + " already imported.", temp_count, "\n -----------")
            
            else:
                print("importing {t_name}".format(t_name=table_name))
                temp_df = pd.read_sql_query(qry, self.conn)
                temp_df.to_csv(self.raw_path+"\\{t_name}.csv".format(t_name=table_name, table=i))
                print("Imported {t_name} succecfully.".format(t_name=table_name), temp_count ,"\n -----------")
            num += 1
            
        print("Done importing database tables.")