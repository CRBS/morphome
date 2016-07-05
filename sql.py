import csv
import sqlite3

def connect(file_sql):
    conn = sqlite3.connect(file_sql)
    cursor = conn.cursor() 
    return conn, cursor

def insert_dataset(cursor, dataset_tuple):
    # Create dataset table if it does not exist
    cursor.execute('CREATE TABLE IF NOT EXISTS dataset '\
        '(id INTEGER PRIMARY KEY, '\
        'ccdb_id INTEGER, '\
        'species TEXT, '\
        'anatomy TEXT, '\
        'zeitgeber_time INTEGER, '\
        'sample_number INTEGER)')

    # Insert new values
    query = 'INSERT INTO dataset (ccdb_id, species, anatomy, zeitgeber_time, '\
        'sample_number) VALUES (?, ?, ?, ?, ?)'
    cursor.execute(query, dataset_tuple)
    return cursor

def insert_cell(cursor, cell_tuple):
    # Create cell table if it does not exist
    cursor.execute('CREATE TABLE IF NOT EXISTS cell '\
        '(id INTEGER PRIMARY KEY, '\
        'dataset_id INTEGER, '\
        'cell_type TEXT, '\
        'FOREIGN KEY(dataset_id) REFERENCES dataset(id))')

    # Insert new values
    query = 'INSERT INTO cell (dataset_id, cell_type) VALUES (?, ?)'
    cursor.execute(query, cell_tuple)
    return cursor

#def insert_organelle(cursor, organelle):
#    # Create organelle table if it does not exist
#    cursor.execute('CREATE TABLE IF NOT EXISTS {0}'\
#        '(id INTEGER PRIMARY KEY, '\
#        'cell_id INTEGER, '\
#        'FOREIGN KEY(cell_id) REFERENCES cell(id))'.format(organelle))
#    return cursor

def get_column_names(cursor, table_name):
    cursor.execute('SELECT * from {0}'.format(table_name))
    column_names = list(map(lambda x: x[0], cursor.description))
    return column_names

def read_csv(cursor, table_name, csv_name, cell_id):
    # Create organelle table if it does not exist
    cursor.execute('CREATE TABLE IF NOT EXISTS {0}'\
        '(id INTEGER PRIMARY KEY, '\
        'cell_id INTEGER, '\
        'FOREIGN KEY(cell_id) REFERENCES cell(id))'.format(table_name))

    existing_column_names = get_column_names(cursor, table_name)
    with open(csv_name, 'r') as fid:
        reader = csv.reader(fid)
        columns = next(reader)
        for col_header in columns:
            if col_header not in existing_column_names:
                query = 'ALTER TABLE {0} ADD {1} FLOAT' 
                query = query.format(table_name, col_header)
                cursor.execute(query)

        query = 'INSERT INTO {0}(cell_id, {1}) VALUES ({2}, {3})'
        query = query.format(table_name, ','.join(columns), cell_id, 
            ','.join('?' * len(columns)))
        for data in reader:
            data = [None if x == 'None' else x for x in data]
            cursor.execute(query, data)
    return cursor  
