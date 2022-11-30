configuration = {
    'pool_size': 10,
    'auto_pcf': {
        'depth': 10000,
        'precision': 50,
        'force_fr': True,
        'timeout_sec': 60,
        'timeout_check_freq': 1000,
        'no_exception': False
    },
    'jobs_to_run': {
        'poly_pslq': {
            'args': { 'degree': (2, 1), 'bulk': 1000, 'subdivide': {
                'PcfCanonical': { 'count': 1, 'balanced_only': True },
                'Named': { 'count': 2 }
                }
            },
            'run_async': True,
            'cooldown': 30,
            'no_work_timeout': 60
        }
    }
}

# If you make your own database, 'name' must match the name in 'create_db.sql' in the line 'CREATE DATABASE <name>'
db_configuration = {
    'host': 'database-1.c1keieal025m.us-east-2.rds.amazonaws.com',
    'port': 5432,
    'user': '',
    'passwd': '',
    'name': 'lirec-main'
}

def get_connection_string(db_name=None):
    conf = db_configuration.copy()
    if db_name:
        conf['name'] = db_name
    return 'postgresql://{user}:{passwd}@{host}:{port}/{name}'.format(**conf)
