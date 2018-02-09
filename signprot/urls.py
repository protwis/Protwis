from django.conf.urls import url

from signprot import views

urlpatterns = [
	url(r'^$', views.BrowseSelection.as_view(), name='index'),
	url(r'^statistics',  views.GProtein, name='gprotein'),
    url(r'^ginterface/(?P<protein>[^/]*?)/$', views.Ginterface, name='render'),
    url(r'^ginterface[/]?$', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^ajax/barcode/(?P<slug>[-\w]+)/(?P<cutoff>\d+\.\d{0,2})/$', views.ajaxBarcode, name='ajaxBarcode'),
    url(r'^ajax/interface/(?P<slug>[-\w]+)/$', views.ajaxInterface, name='ajaxInterface'),
    url(r'^(?P<slug>[-\w]+)/$', views.signprotdetail, name='signprotdetail'),
    url(r'^structure/(?P<pdbname>[-\w]+)/$', views.StructureInfo, name='StructureInfo'),
    url(r'^family/(?P<slug>[-\w]+)/$', views.familyDetail, name='familyDetail'),
    url(r'^matrix$', views.InteractionMatrix, name='InteractionMatrix'),
]
