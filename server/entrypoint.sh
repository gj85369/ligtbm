#!/bin/sh

# wait for the database
echo "Waiting for postgres..."
while ! nc -z ${LOCAL_DB_HOST} 5432; do
  sleep 0.1
done
echo "PostgreSQL started"

# Add changes to the database
python manage.py makemigrations --noinput
python manage.py migrate --noinput
python manage.py migrate django_celery_results --noinput  # celery tutorial says "celery_results", but it gives error

# Copy all the static file to one directory (specified in settings.py)
python manage.py collectstatic --no-input

# Create superuser and an anonymous user
python manage.py shell -c '''
from django.contrib.auth.models import User

if not User.objects.filter(username="admin").exists(): 
	User.objects.create_superuser(username="admin", password="admin", email="")

if not User.objects.filter(username="anonym").exists():
	User.objects.create_user(username="anonym", password="97531anonymous13579", email="")
'''

mkdir -p ${LOCAL_STORAGE}/jobs && mkdir -p ${LOCAL_STORAGE}/tmp

exec "$@"
