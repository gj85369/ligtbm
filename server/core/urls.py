from django.urls import path, include, re_path
from rest_framework import routers
from . import views


router = routers.DefaultRouter()
router.register(r'api/jobs', views.JobViewSet)

urlpatterns = [
    path('', views.index),
    path('index', views.index, name='index'),
    path('login', views.login_page, name='login'),
    path('queue', views.queue_page, name='queue'),
    path('results', views.results_page, name='results'),
    path('publications', views.publications_page, name='publications'),
    path('contact', views.contact_page, name='contact'),
    path('signup', views.signup_page, name='signup'),
    path('logout', views.logout_page, name='logout'),
    path('details', views.details_page, name='details'),
    path('download_file', views.download_file, name='download_file'),
    path('restart_job', views.restart_job, name='restart_job'),
    path('cancel_job', views.cancel_job, name='cancel_job'),
    path('reset_password', views.reset_password_page, name='reset_password'),
    path('settings', views.settings_page, name='settings'),
    path('reject', views.reject_access, name='reject'),
    re_path('^flower/.*', views.flower, name='flower'),
    path('', include(router.urls)),
    path('api/submit/', views.api_submit),
    path('api/download/', views.api_download),
    path('api-auth/', include('rest_framework.urls', namespace='rest_framework')),
]
