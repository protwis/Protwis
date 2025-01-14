"""
Django settings for protwis project.
"""
# Import local settings
# by default, local settings are in protwis/settings_local_development.py
# you can override these settings by creating a protwis/settings_local.py file (or copying settings_local_development)
# protwis/settings_local.py is ignored using .gitignore
try:
    from protwis.settings_local import *
except ImportError:
    from protwis.settings_local_development import *

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))


# Application definition
INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'debug_toolbar',
    'rest_framework',
    'rest_framework_swagger',
    'polymorphic',
    'common',
    'api',
    'news',
    'pages',
    'home',
    'protein',
    'family',
    'residue',
    'alignment',
    'similaritysearch',
    'similaritymatrix',
    'structure',
    'ligand',
    'interaction',
    'mapper',
    'mutation',
    'phylogenetic_trees',
    'sitesearch',
    'build',
    'construct',
    'tools',
    'drugs',
    'signprot',
    'mutational_landscape',
    'contactnetwork',
    'seqsign',
    'angles',
    'hotspots',
    'django_rdkit',
)

MIDDLEWARE = (
    'common.middleware.stats.StatsMiddleware',
    'debug_toolbar.middleware.DebugToolbarMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    # 'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'protwis.urls'

# WSGI_APPLICATION = 'protwis.wsgi.application'

# Max size file upload
FILE_UPLOAD_MAX_MEMORY_SIZE = 26214400

# Internationalization
LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'Europe/Copenhagen'
USE_I18N = True
USE_L10N = True
USE_TZ = True

# Default site configuration (gpcr - GPCRdb, gprotein - GproteinDb, arrestin - ArrestinDb, biasedsignalingatlas - Biased Signaling Atlas)
DEFAULT_SITE = "gpcr"

# Static files (CSS, JavaScript, Images)
STATIC_URL = '/static/'
STATIC_ROOT = '/protwis/static/protwis'
STATICFILES_DIRS = (os.sep.join([BASE_DIR, "static"]),)
MEDIA_URL = '/media/'
MEDIA_ROOT = '/protwis/media/protwis'

# Serializer
SESSION_SERIALIZER = 'django.contrib.sessions.serializers.PickleSerializer'
SESSION_COOKIE_AGE = 86400 #Expire cookies and session after 24 hrs
SWAGGER_SETTINGS = {
    'USE_SESSION_AUTH' : False,
    "exclude_namespaces": ["excluded_apis"],
}

# Templates
TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            # insert your TEMPLATE_DIRS here
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',
                'protwis.context_processors.current_site',
                'protwis.context_processors.canonical_tag',
                'protwis.context_processors.documentation_url',
                'protwis.context_processors.google_analytics',
                'protwis.context_processors.site_title'
            ],
        },
    },
]

if DEBUG:
    TEMPLATES[0]['OPTIONS']['debug'] = True

# Debug toolbar
if DEBUG:
    DEBUG_TOOLBAR_PATCH_SETTINGS = False
    INTERNAL_IPS = ('10.0.2.2')

SMILES_MAX_LENGTH = 200
PYTHON_SMILES_VALIDATION = True

# Logging
if DEBUG:
    LOGGING = {
       'version': 1,
       'disable_existing_loggers': False,
       'formatters': {
           'verbose': {
               'format' : "[%(asctime)s] %(levelname)s [%(name)s:%(lineno)s] %(message)s",
               'datefmt' : "%d/%b/%Y %H:%M:%S"
           },
       },
       'handlers': {
           'django': {
               'level': 'WARNING',
               'class': 'logging.FileHandler',
               'filename': 'logs/django.log',
               'formatter': 'verbose'
           },
           'build': {
               'level': 'DEBUG',
               'class': 'logging.FileHandler',
               'filename': 'logs/build.log',
               'formatter': 'verbose'
           },
           'protwis': {
               'level': 'DEBUG',
               'class': 'logging.FileHandler',
               'filename': 'logs/protwis.log',
               'formatter': 'verbose'
           },
       },
       'loggers': {
           'django': {
               'handlers':['django'],
               'propagate': True,
               'level':'WARNING',
           },
           'build': {
               'handlers': ['build'],
               'level': 'DEBUG',
           },
           'protwis': {
               'handlers': ['protwis'],
               'level': 'DEBUG',
           },
       }
    }

#CACHE
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
        'LOCATION': '/tmp/django_cache',
        'OPTIONS': {
            'MAX_ENTRIES': 10000000
        }
    },
    'alignments': {
        'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
        'LOCATION': '/tmp/django_cache_alignments',
        'OPTIONS': {
            'MAX_ENTRIES': 1000
        }
    }
}

# Note that https://www.django-rest-framework.org/community/3.10-announcement
# So, have to switch from CoreAPI to OpenAPI. Next line will work for now.
# Uncomment when needed.
# REST_FRAMEWORK = { 'DEFAULT_SCHEMA_CLASS': 'rest_framework.schemas.coreapi.AutoSchema' }
