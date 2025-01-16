from django.conf import settings
from django.core.cache import cache

from common.models import Citation

def current_site(request):
    domain_switches = {"gpcrdb.org" : "gpcr", "gproteindb.org" : "gprotein", "arrestindb.org": "arrestin", "biasedsignalingatlas.org": "biasedsignalingatlas"}
    domain = request.get_host().lower()

    if not domain in domain_switches:
        return {
           'current_site': settings.DEFAULT_SITE
         }
    else:
        return {
           'current_site': domain_switches[domain]
         }

def site_title(request):
    domain = current_site(request)["current_site"]
    domain_titles = {"gpcr": "GPCRdb", "gprotein": "GproteinDb", "arrestin": "ArrestinDb", "biasedsignalingatlas": "Biased Signaling Atlas"}

    if not domain in domain_titles:
        domain = settings.DEFAULT_SITE

    return {
       'site_title': domain_titles[domain]
     }

def canonical_tag(request):
    citation_dict = cache.get("citation_urls")
    if citation_dict == None:
        citation_dict = {}
        citation_urls = Citation.objects.all().values_list("url", flat = True)
        for url in citation_urls:
            path = url.split(".org")[1]
            citation_dict[path] = url
        cache.set("citation_urls", citation_dict, 60*60*24*7)

    if request.path in citation_dict:
        return {
           'canonical_tag': citation_dict[request.path]
         }
    elif request.path == "" or request.path == "/":
        return {
           'canonical_tag': "https://" + request.get_host()
         }
    else:
        return {
           'canonical_tag': "https://gpcrdb.org" + request.path
         }

def documentation_url(request):
    return {
        'documentation_url': settings.DOCUMENTATION_URL
    }

def google_analytics(request):
    """
    Use the variables returned in this function to
    render your Google Analytics tracking code template.
    """
    if settings.GOOGLE_ANALYTICS_KEY:
        return {
            'google_analytics': settings.GOOGLE_ANALYTICS_KEY
        }
    return {}
