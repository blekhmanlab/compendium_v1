# Params for identifying ourselves to the NCBI API
tool = 'name_for_auditing'
email = 'your_contact_here'

# DB connection info. This should be simple because it's local
db = {
    'host': 'localhost',
    'db': 'compendiumdb',
    'user': 'postgres',
    'password': '',
    'schema': 'full',
    'connection': {
        'max_attempts': 5,
        'timeout': 20,
        'attempt_pause': 5
    }
}
