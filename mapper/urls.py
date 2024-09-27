from django.conf.urls import url
from django.views.decorators.cache import cache_page
from . import views

urlpatterns = [
    url(r'^$', views.LandingPage.as_view(), name='landing_page'),
    url(r'^plotrender', views.plotrender.as_view(), name='data_mapper_plotrender')
]
