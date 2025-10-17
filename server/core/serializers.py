import logging
from rest_framework import serializers
from .models import Job

logger = logging.getLogger('core')


class JobSerializer(serializers.ModelSerializer):
    class Meta:
        model = Job
        fields = ('job_id', 'job_name', 'status')

    def update(self, instance, validated_data):
        raise RuntimeError('Not implemented')
