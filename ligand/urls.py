from django.conf.urls import url
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from django.db import models
from ligand.views import *

urlpatterns = [
    url(r'^$',cache_page(3600*24*7)(LigandBrowser.as_view()), name='ligand_browser'),
    url(r'^target/compact/(?P<pk>[-\w]+)$',TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^target/(?P<pk>[-\w]+)$',TargetDetails, name='ligand_target_detail'),
    url(r'^target/purchasable/(?P<pk>[-\w]+)$',TargetPurchasabilityDetails, name='ligand_target_detail_purchasable'),
    url(r'^(?P<ligand_id>[-\w]+)/$',LigandDetails, name='ligand_detail'),
    url(r'^statistics', cache_page(3600*24*7)(LigandStatistics.as_view()), name='ligand_statistics'),
    url(r'^experiment/(?P<pk>[-\w]+)/detail$', ExperimentEntryView.as_view()),
    url(r'^vendors$', test_link, name='test'),
    # url(r'^biasedbrowser$', cache_page(3600*24*7)(BiasBrowser.as_view()), name='bias_browser'),
    # url(r'^biasedgbrowser$', cache_page(3600*24*7)(BiasBrowserGSubbtype.as_view()), name='bias_g_browser'),
    url(r'^biasedbrowser$', BiasBrowser.as_view(), name='bias_browser'),
    url(r'^biasedgbrowser$', BiasBrowserGSubbtype.as_view(), name='bias_g_browser'),
    url(r'^browserchembl$', BiasBrowserChembl.as_view(), name='bias_chembl_browser'),
    url(r'^browservendors$', BiasVendorBrowser.as_view(), name='browservendor'),
    url(r'^biasedpathways$', BiasPathways.as_view(), name='pathways'),
    url(r'^pathwaydata/(?P<pk>[-\w]+)/detail$', PathwayExperimentEntryView.as_view()),
    url(r'^(?P<pk>[-\w]+)/info$', LigandInformationView.as_view()),
    # url(r'^targetselection', BiasTargetSelection.as_view(), name='targetselection'),
    # url(r'^gtargetselection', BiasGTargetSelection.as_view(), name='gtargetselection'),

]
