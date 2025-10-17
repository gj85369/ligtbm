import os

CONFIG = {
    "database": {
        'NAME': os.getenv('LOCAL_DB_NAME', 'server_template'),
        'HOST': os.getenv('LOCAL_DB_HOST', '127.0.0.1'),
        'USER': os.getenv('LOCAL_DB_USER', 'postgres'),
        'PASS': None,
        'PORT': 5432
    },
    "email": {
        'HOST': os.getenv('EMAIL_HOST'),
        'USER': os.getenv('EMAIL_USER'),
        'PASS': os.getenv('EMAIL_PASS'),
        'PORT': 587
    },
    "remote": {
        'HOST': os.getenv('REMOTE_HOST'),
        'USER': os.getenv('REMOTE_USER'),
        'PASS': os.getenv('REMOTE_PASS'),
        'STORAGE': os.getenv('REMOTE_STORAGE'),
        'BIN': os.getenv('REMOTE_BIN'),
    },
    "local": {
        'ROOT': os.getenv('LOCAL_ROOT'),
        'STORAGE': os.getenv('LOCAL_STORAGE'),
        'SECRET_KEY': os.getenv('SECRET_KEY')
    },
    "rabbit": {
        'USER': os.getenv('RABBITMQ_DEFAULT_USER'),
        'PASS': os.getenv('RABBITMQ_DEFAULT_PASS')
    },
    "production": (os.getenv('PRODUCTION', '0').lower() in ['true', 'yes', '1'])
}
