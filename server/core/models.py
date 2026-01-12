from django.db import models
from django.contrib.auth.models import User
from django.core.files.base import File
from . import env
from .utils import upload_file
from .runner.runner import run_job

from django.core.validators import (
    FileExtensionValidator,
    MinLengthValidator,
    RegexValidator,
)
from .modules.validators import validate_scc_api_url

from celery import shared_task
#from celery.utils.log import get_task_logger
from django.db.models import (
    CASCADE,
    BooleanField,
    CharField,
    DateTimeField,
    FileField,
    FloatField,
    ForeignKey,
    IntegerChoices,
    IntegerField,
    Model,
    TextChoices,
    TextField,
    UUIDField,
)
import json
import logging
import tempfile
from time import sleep
from path import Path


def _get_task_logger(job_id):
    logger = logging.getLogger(str(job_id))
    logger.setLevel(logging.DEBUG)
    job = Job.objects.get(job_id=job_id)
    h = logging.FileHandler(job.get_dir().joinpath('task.log'), mode='a')
    f = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s')
    h.setFormatter(f)
    logger.addHandler(h)
    return logger


@shared_task(bind=True, acks_late=True, track_started=True, default_retry_delay=30, autoretry_for=(Exception,), max_retries=1)
def submit_job(self, job_id):
    logger = _get_task_logger(job_id)
    try:
        job = Job.objects.get(job_id=job_id)
        cur_status = job.status
        if cur_status in Job.STATUS_FINISHED:
            logger.info(f'Restarting job {job_id} in L.STR (currently in {cur_status})')
            job.status = 'L.STR'
            job.save()
        else:
            logger.info(f'Starting job {job_id} in {cur_status} (it was interrupted half-way)')
        run_job(job)
    except Exception as e:
        logger.exception(e)
        Job.objects.get(job_id=job_id).local_error('Unknown error')
        raise

class Settings(Model):                                                                                                                                                                        
    ACCEPTING_JOBS = BooleanField(default=False, help_text="Accepting jobs if True.")                                                                                                         
    BANNER_ON = BooleanField(default=False, help_text="Banner on if True.")                                                                                                                   
    BANNER_MESSAGE = TextField(                                                                                                                                                               
        blank=True,                                                                                                                                                                           
        null=True,                                                                                                                                                                            
        help_text="Enter a message to be displayed on the site-wide banner.",                                                                                                                 
    )                                                                                                                                                                                         
    AUTO_BANNER_MESSAGE = TextField(                                                                                                                                                          
        default="Submissions are temporarily shut off due to connection errors with the Shared Computing Center. We are working to resolve this as quickly as possible.",                     
        help_text="This is the text that will be automatically displayed when the connection with the Vajda-Dashboard is interrupted.",                                                       
    )                                                                                                                                                                                         
    DELETE_JOBS_AFTER = IntegerField(                                                                                                                                                         
        default=30, help_text="Delete jobs that have not been modified in X days."                                                                                                            
    )                                                                                                                                                                                         
    TIMEOUT_JOBS_AFTER = IntegerField(                                                                                                                                                        
        default=30, help_text="Timeout jobs that have not been modified in X hours."                                                                                                          
    )                                                                                                                                                                                         
    SCC_API_URL = CharField(                                                                                                                                                                  
        max_length=50,                                                                                                                                                                        
        validators=[validate_scc_api_url],                                                                                                                                                    
        help_text="Enter the vajda-dashboard API URL",                                                                                                                                        
    )                                                                                                                                                                                         
    SCC_API_TOKEN = CharField(                                                                                                                                                                
        primary_key=True, max_length=40, validators=[MinLengthValidator(40)]                                                                                                                  
    )                                                                                                                                                                                         
    modified = DateTimeField(auto_now=True)                                                                                                                                                   

class SettingsLog(Model):                                                                                                                                                                     
    settings = ForeignKey("Settings", on_delete=CASCADE)                                                                                                                                      
    event = TextField()                                                                                                                                                                       
    created = DateTimeField(auto_now_add=True)                                                                                                                                                
                                                                                                                                                                                              
    class Meta:                                                                                                                                                                               
        get_latest_by = ["created"]                                                                                                                                                           
        ordering = ["-created"]                                                                                                                                                               
                                                                                                                                                                                              
    def __str__(self):                                                                                                                                                                        
        """ return string with the log pk and the associated event"""                                                                                                                         
        return f"{self.pk}: {self.event}"    

        
class Job(models.Model):
    STATUS_CHOICES = (
        ('L.STR', 'Starting the job'),
        ('R.RUN', 'Running'),
        ('L.FIN', 'Finalizing job'),
        ('L.CPL', 'Job completed'),
        ('L.PDB', 'No templates found'),
        ('L.ERR', 'Error'),
        ('R.ERR', 'Error'),
    )

    # All other statuses mean the job is still running
    STATUS_FINISHED = ('L.CPL', 'L.PDB', 'L.ERR', 'R.ERR')

    # SCC queue status
    QUEUE_CHOICES = (
        ('QUE', 'in queue'),
        ('RUN', 'running'),
        ('ERR', 'exited with error'),
        ('FIN', 'finished'),
        ('NAN', '')
    )

    job_id = models.AutoField(primary_key=True)
    job_name = models.CharField(max_length=100, default="", blank=True)
    created = models.DateTimeField(auto_now_add=True)
    warning = models.CharField(max_length=1000, default="", blank=True)
    error = models.CharField(max_length=1000, default="", blank=True)
    deleted = models.BooleanField(default=False)
    restarted = models.IntegerField(default=0)
    touched = models.DateTimeField(auto_now=True)
    queue_id = models.IntegerField(blank=True, null=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    ip = models.GenericIPAddressField()
    celery_id = models.CharField(max_length=200, default="", blank=True)

    status = models.CharField(max_length=5, choices=STATUS_CHOICES, default='L.STR')
    queue_status = models.CharField(max_length=3, choices=QUEUE_CHOICES, default='NAN')
    details_json = models.CharField(max_length=1000, blank=True)

    _temp_dir = ''

    def __str__(self):
        attr_list = [self.job_id, self.user, self.status]
        return '-'.join(map(str, attr_list))

    def __getattr__(self, attr):
        json_str = super(Job, self).__getattribute__("details_json")
        if json_str:
            json_dict = json.loads(json_str)
            if attr in json_dict:
                return json_dict[attr]
        return super(Job, self).__getattribute__(attr)

    def get_dir(self):
        return env.JOBS_DIR.joinpath(str(self.job_id))

    def get_remote_dir(self):
        return env.REMOTE_JOBS.joinpath(str(self.job_id))

    def get_output_dir(self):
        return self.get_dir().joinpath('output')

    def get_user_dir(self):
        return self.get_dir().joinpath('user')

    def get_user_json(self):
        return self.get_user_dir().joinpath('user.json')

    def get_models_dir(self):
        return self.get_dir().joinpath('minimized')

    def get_scores_csv(self):
        return self.get_output_dir().joinpath('scores.csv')

    def get_models_zip(self):
        return self.get_output_dir().joinpath('models.zip')

    def get_log_txt(self):
        return self.get_output_dir().joinpath('log.txt')

    def get_receptor_pdb(self):
        return self.get_user_dir().joinpath(self.rec_pdb_file)

    def remote_error(self, error):
        self.status = 'R.ERR'
        self.error = error
        self.save()

    def local_error(self, error):
        self.status = 'L.ERR'
        self.error = error
        self.save()

    def reset(self):
        """
        Call this when a job is (re)started. Resets error and warning and increments 'restarted' field
        """
        self.restarted += 1
        self.error = ''
        self.warning = ''
        self.save()

    def start(self):
        if self.celery_id != '':
            result = submit_job.AsyncResult(self.celery_id)
            if not result.ready():
                return [f'Job {self.job_id} is already running (state: {result.state})']
        result = submit_job.apply_async((self.job_id,))
        self.celery_id = result.task_id
        self.restarted += 1
        self.error = ''
        self.warning = ''
        self.save()
        return []

    def cancel(self):
        if self.celery_id != '':
            result = submit_job.AsyncResult(self.celery_id)
            if result.ready():
                return [f'Job has already finished (state: {result.state})']

            result.revoke(terminate=True)
            return [f'Job cancelled (state: {result.state})']
        else:
            return ['Task id is empty']

    def _make_temp_dir(self):
        self._temp_dir = Path(tempfile.mkdtemp(dir=env.TMP_DIR))
        self._temp_dir.chmod(0o755)

    def _create_job_dir(self):
        job_dir = self.get_dir()
        if job_dir.isdir():
            job_dir.rmtree(job_dir)
        self._temp_dir.move(job_dir)

    @staticmethod
    def _get_name_for_uploaded_file(field_name, file):
        src_name = file.name
        src_ext = Path(src_name).ext  # Python's standard `pathlib` uses .suffix, but `path` uses .ext
        out_name = Path(field_name).with_suffix(src_ext)
        return out_name

    def check_user_input(self, form, files):
        error_list = []
        self._make_temp_dir()

        user_dir = self._temp_dir / 'user'
        user_dir.mkdir_p()

        # dump user files in user/
        for name, file in files.items():
            fpath = user_dir / self._get_name_for_uploaded_file(name, file)
            upload_file(file, fpath)

        # dump user data in user/user.json
        with open(user_dir.joinpath('user.json'), 'w') as f:
            data = {k: (
                v if not isinstance(v, File)
                else self._get_name_for_uploaded_file(k, v)) for k, v in form.cleaned_data.items()}
            json.dump(data, f, indent=4)

        return error_list

    def make_dir_and_start(self):
        self._create_job_dir()

        clean_json = self.get_user_json()
        with open(clean_json, 'r') as f:
            details_json = json.load(f)
        self.details_json = json.dumps(details_json)
        self.save()

        self.start()


class ModellerKey(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE, primary_key=True)
    key = models.CharField(max_length=300, blank=True, null=True)
