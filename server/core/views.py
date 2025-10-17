import logging

from django.http import HttpResponse, Http404
from django.shortcuts import render, redirect, reverse
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required, user_passes_test
from django.views.decorators.csrf import csrf_exempt
from django.core.paginator import Paginator
from django.core.exceptions import PermissionDenied
from django.db import IntegrityError
from django.contrib.auth.hashers import make_password
from rest_framework import viewsets, permissions, status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from .models import Job, ModellerKey
from .serializers import JobSerializer
from . import env
from . import forms
from . import emails
from . import utils

logger = logging.getLogger('core')


def _get_client_ip(request):
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
    if x_forwarded_for:
        ip = x_forwarded_for.split(',')[0]
    else:
        ip = request.META.get('REMOTE_ADDR')
    return ip


def _get_modeller_key(request):
    modeller_key = ''
    if request.user.username != 'anonym':
        try:
            key = ModellerKey.objects.get(user=request.user)
            modeller_key = key.key
        except:
            pass
    return modeller_key


@login_required(login_url='login')
def index(request):
    form = forms.JobSubmitForm()
    errors = []

    if request.method == 'POST':
        form = forms.JobSubmitForm(request.POST, request.FILES)

        if not form.is_valid():
            err_dict = form.errors.get_json_data(escape_html=False)
            errors += [str(v[0]['message']) for k, v in err_dict.items()]
        else:
            # This stage we need, if there is some more complicated processing,
            # for which a temporary directory is required. So a temp directory
            # is created in storage/tmp and you can add some extra preparation/validation
            # there. If there is no errors this temp dir is moved to storage/jobs and
            # becomes permanent, otherwise it gets erased and user receives an error.
            job = Job(user=request.user, ip=_get_client_ip(request), job_name=form.cleaned_data['job_name'])
            errors += job.check_user_input(form, request.FILES)

            if not errors:
                # save modeller key if provided
                if form.cleaned_data['modeller_key'] != '':
                    ModellerKey(user=request.user, key=form.cleaned_data['modeller_key']).save()

                # start job
                job.save()
                job.make_dir_and_start()
                return redirect('queue')

    context = {'errors': errors,
               'form': form,
               'modeller_key': _get_modeller_key(request)}

    return render(request, 'core/home.html', context)


@user_passes_test(lambda u: u.is_superuser, 'reject')
@login_required(login_url='login')
def restart_job(request):
    job_id = int(request.GET.get('job_id'))
    job = Job.objects.get(job_id=job_id)
    errors = job.start()

    if not errors:
        return redirect('queue')

    return HttpResponse('<h1>' + errors[0] + '</h1>')


@user_passes_test(lambda u: u.is_superuser, 'reject')
@login_required(login_url='login')
def cancel_job(request):
    job_id = int(request.GET.get('job_id'))
    job = Job.objects.get(job_id=job_id)
    errors = job.cancel()

    if not errors:
        return redirect('queue')

    return HttpResponse('<h1>' + errors[0] + '</h1>')


def _jobs_page(request, job_list):
    if not request.user.is_superuser:
        job_list = job_list.filter(user__username=request.user.username)

    order_by = request.GET.get('order_by', '-job_id')
    paginator = Paginator(job_list.order_by(order_by, '-job_id'), env.JOBS_PER_PAGE)
    page = request.GET.get('page')
    jobs = paginator.get_page(page)
    return render(request, 'core/jobs.html', {'jobs': jobs, 'order_by': order_by})


@login_required(login_url='login')
def queue_page(request):
    job_list = Job.objects.exclude(status__in=Job.STATUS_FINISHED)
    return _jobs_page(request, job_list)


@login_required(login_url='login')
def results_page(request):
    job_list = Job.objects.filter(status__in=Job.STATUS_FINISHED)
    return _jobs_page(request, job_list)


@login_required(login_url='login')
def details_page(request):
    job_id = request.GET.get('job_id')
    try:
        job = Job.objects.get(job_id=job_id)
    except Job.DoesNotExist:
        reject_access(request)
        return
    except ValueError:
        return HttpResponse(status=status.HTTP_422_UNPROCESSABLE_ENTITY, content='Invalid job_id: should be integer')

    user = request.user

    if not request.user.is_superuser and job.user.id != user.id:
        reject_access(request)

    return render(request, 'core/details.html', {'job': job})


def _make_response_from_bytes(file_data, presented_file_name, content_type='plain/text'):
    resp = HttpResponse(content_type=content_type)
    resp['Content-Disposition'] = 'attachment; filename="{}"'.format(presented_file_name)
    resp.write(file_data)
    return resp


def _make_response_from_file(file_name, presented_file_name, content_type='plain/text'):
    with open(file_name, 'rb') as f:
        return _make_response_from_bytes(f.read(), presented_file_name, content_type)


@login_required(login_url='login')
def download_file(request):
    user = request.user
    job_id = request.GET.get('job_id')
    content = request.GET.get('content')
    try:
        job = Job.objects.get(job_id=job_id)
    except Job.DoesNotExist:
        raise Http404('Wrong Job ID')
    except ValueError:
        return HttpResponse(status=status.HTTP_422_UNPROCESSABLE_ENTITY, content='Invalid job_id: should be integer')

    if not request.user.is_superuser and job.user.id != user.id:
        reject_access(request)

    try:
        if content == 'ngl_models.pdb':
            response = _make_response_from_bytes(utils.concat_models_zipped(job.get_models_zip()),
                                                 content)
        elif content == 'models.zip':
            response = _make_response_from_file(job.get_models_zip(), content, 'application/zip')
        elif content == 'scores.csv':
            response = _make_response_from_file(job.get_scores_csv(), content)
        elif content == 'log.txt':
            response = _make_response_from_file(job.get_log_txt(), content)
        elif content == 'receptor.pdb':
            response = _make_response_from_file(job.get_receptor_pdb(), content)
        else:
            raise Http404('Wrong content type')
    except FileNotFoundError:
        raise Http404('Content does not exist for this job')

    return response


def reject_access(request):
    raise PermissionDenied('You don\'t have permission to access this page')


def publications_page(request):
    return render(request, 'core/publications.html')


def contact_page(request):
    return render(request, 'core/contact.html')


def signup_page(request):  # TODO: add second password field
    logout(request)

    context = {'page_title': 'Sign Up',
               'page_description': '(you must provide your academic e-mail address)',
               'submit_text': 'Sign Up'}

    errors = []

    if request.method == 'POST':
        form = forms.SignUpForm(request.POST)

        if not form.is_valid():
            err_dict = form.errors.get_json_data(escape_html=False)
            errors += [str(v[0]['message']) for k, v in err_dict.items()]
        else:
            new_password = utils.random_string(10)
            user = User(username=form.cleaned_data['email'],
                        password=make_password(new_password),
                        first_name=form.cleaned_data['first_name'],
                        last_name=form.cleaned_data['last_name'],
                        email=form.cleaned_data['email'])
            try:
                user.save()
            except IntegrityError as e:
                errors.append(f'Failed to create a user, please <a href={reverse("contact")}>contact us</a>')
            else:
                # you can get SMTPAuthenticationError if you didnt switch your google account to "less secure apps" mode
                emails.send_greeting(user, new_password)
                context.update({'email': user.email})
                return render(request, 'core/thankyou.html', context)

    else:
        form = forms.SignUpForm()

    context.update({'form': form, 'errors': errors})
    return render(request, 'core/generic_form.html', context)


def login_page(request):
    logout(request)

    if request.GET.get('anonym') == 'true':
        username = 'anonym'
        password = '97531anonymous13579'
        user = authenticate(request, username=username, password=password)

        if user is not None:
            login(request, user)
            return redirect('index')
        else:
            return render(request, 'core/login.html', {'errors': ['Cannot use anonymous login due to internal error']})

    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)

        if user is not None:
            login(request, user)
            return redirect('index')
        else:
            return render(request, 'core/login.html', {'errors': ['Invalid username or password']})

    return render(request, 'core/login.html')


def logout_page(request):
    logout(request)
    return redirect('login')


def reset_password_page(request):
    logout(request)

    context = {'page_title': 'Reset password',
               'page_description': 'We will send your new password to the provided e-mail address',
               'submit_text': 'Reset password'}

    errors = []
    messages = []
    form = forms.PasswordResetForm()

    if request.method == 'POST':
        form = forms.PasswordResetForm(request.POST)

        if not form.is_valid():
            err_dict = form.errors.get_json_data(escape_html=False)
            errors += ['{}'.format(v[0]['message']) for k, v in err_dict.items()]
        else:
            try:
                new_password = utils.random_string(10)
                username = form.cleaned_data['username']
                user = User.objects.get(username=username)
                user.set_password(new_password)
                emails.send_password(user, new_password)
            except Exception as e:
                logging.exception(e)
                errors += [f'Internal error has occurred. Please <a href={reverse("contact")}>contact us</a> to solve the issue.']
            else:
                user.save()
                messages += ['Your password was successfully reset, please check your e-mail.']
                return render(request, 'core/login.html', {'messages': messages, 'form': form})

    context.update({'form': form, 'errors': errors})
    return render(request, 'core/generic_form.html', context)


@user_passes_test(lambda u: u.username != 'anonym', 'reject')
@login_required(login_url='login')
def settings_page(request):
    context = {'page_title': 'Settings',
               'page_description': 'Change your account preferences here',
               'submit_text': 'Save changes'}

    errors = []
    messages = []
    form = forms.SettingsForm()
    if request.method == 'POST':
        form = forms.SettingsForm(request.POST)
        if not form.is_valid():
            err_dict = form.errors.get_json_data(escape_html=False)
            errors += [str(v[0]['message']) for k, v in err_dict.items()]
        else:
            if not request.user.check_password(form.cleaned_data['current_password']):
                errors += ['Wrong password provided']
                form.add_error('current_password', 'Wrong password provided')
            else:
                new_password = form.cleaned_data['password']
                try:
                    request.user.set_password(new_password)
                    request.user.save()
                except Exception as e:
                    logger.exception(e)
                    errors += ['Internal error occurred']
                else:
                    messages += ['Your changes were successfully applied']
                    context.update({'form': form, 'messages': messages})
                    return render(request, 'core/generic_form.html', context)

    context.update({'form': form, 'errors': errors})
    return render(request, 'core/generic_form.html', context)


@csrf_exempt  # TODO: Let's hope this is safe. Otherwise django doesn't let flower make ajax requests
@user_passes_test(lambda u: u.is_superuser, 'reject')
@login_required(login_url='login')
def flower(request):
    response = HttpResponse()
    path = request.get_full_path()
    path = path.replace('flower', '_flower_internal', 1)
    response['X-Accel-Redirect'] = path
    return response


class JobViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Job.objects.all()
    serializer_class = JobSerializer
    permission_classes = [permissions.IsAdminUser]


@api_view(['POST'])
def api_submit(request):
    permission_classes = [permissions.IsAdminUser]
    form = forms.JobSubmitForm(request.POST, request.FILES)
    if not form.is_valid():
        return Response({'status': 'L.ERR', 'errors': form.errors.get_json_data(escape_html=True)}, status=status.HTTP_400_BAD_REQUEST)

    job = Job(user=request.user, ip=_get_client_ip(request), job_name=form.cleaned_data['job_name'])
    errs = job.check_user_input(form, request.FILES)

    if errs:
        return Response({'status': 'L.ERR', 'errors': form.errors.get_json_data(escape_html=True)}, status=status.HTTP_400_BAD_REQUEST)
    job.save()

    # save modeller key if provided
    if form.cleaned_data['modeller_key'] != '':
        ModellerKey(user=request.user, key=form.cleaned_data['modeller_key']).save()

    # start job
    job.make_dir_and_start()

    return Response({'job_id': job.job_id, 'status': job.status}, status=status.HTTP_201_CREATED)


@api_view(['GET'])
def api_download(request):
    permission_classes = [permissions.IsAdminUser]
    response = download_file(request)
    if response.status_code == 200:
        return response
    else:
        return Response({'errors': response.content},
                        status=response.status_code)
