# LigTBM
#### A server for template-based small molecule docking
Powered by Django + Nginx + Postgres + Bootstrap + Celery.

------
### Getting started
1. [Install Docker Community Edition](https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/)
2. Install docker-compose into python3, e.g. `pip3 install --user docker-compose`
3. Add your user to the docker group. `sudo usermod -a -G docker username`;
you may have to reboot after this step for you to show up in the group.

You should then use the `local-docker-compose` script as a drop in replacement
for docker-compose. For example, to start the server you can 
run `local-docker-compose up --build`. Before you do that you need to create `.local_params`
in the root directory. Read below on how to do this.

Cleaning up after docker for a clean rebuild:    
1. `./cluspro-docker-compose rm` will remove the containers    
2. `docker volume prune`

If you don't explicitly remove the volumes between docker runs, the databases persist, 
so you can stop the containers and launch them again safely without any loss of data.

Before you can run the jobs, you need to [deploy the backend on SCC](https://bitbucket.org/abc-group/ligtbm/src/master/backend)

#### Architecture
Docker runs several services: web (which runs Gunicorn), nginx, db (Postgres database), rabbitmq, 
celery, and flower. 
Gunicorn handles the python (Django) code, accesses the database and cooperates with Nginx.
Celery is a background task manager and it needs RabbitMQ to run (message broker). Flower is
a task monitor, which is powered by Celery.

The production configuration is described in greater details in [a separate document](https://docs.google.com/document/d/1gSqF5NH0GEjX5KystSOK1FhP4nK_PRGnpZrRcda9rCg/).

#### Structure
All the frontend code is located in `server/`.
The structure of `server/` directory is enforced by the Django rules, so we have
`server/server`, where all the server settings are located (`settings.py`) as well as `config.py`. 
`config.py` is where the custom variables are kept (e.g. email login 
 and password for sending messages to the user), which in turn are populated from 
 the environment, which is set in docker-compose.yml.     
`core/` contains the app code, as it's called in Django. `core/templates` has all 
the html files, `core/static` - CSS and JS, and `runner/` contains the code for job 
running.

#### `Core/` structure
1. `views.py` is the main file - it contains all the page renderers and handles 
all the forms and requests. Most of functions return an HTML response.

2. `urls.py` assigns URLs to the functions in `views.py`.

3. `models.py` contains custom data tables, which are added to the default Django
tables. Right now it contains a model for jobs, which can be customized as you wish.
The intention, however, was to keep all the generic job fields as separate class attributes
(job name, IP etc.) and to store all the rest job-type specific parameters as a json string
in details_json field. This way we can prevent creating many different tables for different 
job types or addition of infinite new fields to the same job table (once we add new job parameters, 
for example).

4. All the forms, which are on the website are contained in `forms.py` and it should
be kept so. These forms are all handled in `views.py`.

5. `emails.py` has messages for users, whenever we want to send them something. They
use the e-mail address and password specified in `server/settings.py`, which are in 
turn taken from environmental variables in `docker-compose.yml`. If they were not 
specified you will get an error, whenever the server is trying to send a message.

6. `env.py` is where you should keep your local variables. Also all the variables
 in `env` dictionary will be passed as a context to the html templates, so you 
 can refer to them in the templates.
 
#### At the first launch
Before the first launch, please create `.local_params` file, using `.local_params_example` as a template.
Some of the parameters are discussed below.
Upon launch, the database is initialized and two users are created:     
1) admin with password 'admin'. This is a superuser, you should change the 
password for it immediately. The admin page is located at http://localhost/admin    
2) anon, which is where you log in once you click 'use without your own account' 
button on the login page. It has limited permissions.    
The `storage/` directory is created to keep the jobs.

#### Jobs
When you run jobs they are stored in docker container in `/storage`, which is 
by default mounted in your project root. You can change this by adding `LOCAL_STORAGE` 
parameter to `.local_params`.
Storage has two directories: `tmp/` for temporary storage, if you need to compute
or check something before adding the job to the database, and `jobs/` with all
the jobs.

-----
### Running jobs
#### Backend
Deploy the backend on SCC before running jobs.

#### .local_params
Environmental variables with some paths, e-mail login and password are stored 
in `.local_params`, which are used when you run `local-docker-compose`. To create the 
file use `.local_params_example` as a template.

`REMOTE_HOST` - SCC server name (default: scc2.bu.edu)    
`REMOTE_USER` - your user on SСС    
`REMOTE_PASS` - your password    
`REMOTE_STORAGE` - where you want to keep the jobs on SCC (must exist)    
`REMOTE_BIN` - where you deployed the backend on SCC    

Variables for sending e-mails. If you don't specify them you will get errors when new users are registered etc.
If your e-mail is `server@gmail.com` and the password is `password` then the values should be:

`EMAIL_USER` - server    
`EMAIL_PASS` - password     
`EMAIL_HOST` - smtp.gmail.com

`LOCAL_PORT` is the port through which you access the server. We recommend to also specify IP address
to listen on, e.g., `127.0.0.1:8080` to listen on local port 8080. If only port number is specified,
the Docker will allow connections from anywhere in the world.

`SECRET_KEY` is for Django internal use (is generated at the first run 
of `local-docker-compose`) and should be kept secret.

`RABBITMQ_USER` is the username for accessing built-in RabbitMQ server. As the server is, by default, 
only accessible from other docker containers within this server, it could safely be kept default.    
`RABBITMQ_PASS` is the password for the built-in RabbitMQ server. If not specified, it will be
 autogenerated.
 
 `LOCAL_STORAGE` is the directory on the filesystem where job data is stored (see "Jobs" section).     
 `DATABASE_VOLUME` is the directory on the filesystem where database is stored. If not specified,
 a dedicated docker container `postgres_data` is used.

`PRODUCTION` - set to `1` to use production settings for `ligtbm.cluspro.org` (hostname filter, no debug output, 
emailing error reports).