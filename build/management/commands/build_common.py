from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from common.models import WebResource, WebLink, PublicationJournal, Publication
from common.tools import test_model_updates
from protein.models import (ProteinSegment, ProteinAnomaly, ProteinAnomalyType, ProteinAnomalyRuleSet,
    ProteinAnomalyRule, Site)
from ligand.models import Ligand, LigandType, LigandRole
from residue.models import ResidueGenericNumber, ResidueNumberingScheme
from news.models import News

import django.apps
import logging
import shlex
import os
import yaml

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)
    resource_source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'resources.txt'])
    ligands_source_file = os.sep.join([settings.DATA_DIR, 'ligand_data', 'ligands.yaml'])
    publications_source_file = os.sep.join([settings.DATA_DIR, 'publications_data', 'publications.yaml'])
    segment_source_file = os.sep.join([settings.DATA_DIR, 'protein_data', 'segments.txt'])
    residue_number_scheme_source_file = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers',
        'schemes.txt'])
    anomaly_source_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'anomalies'])
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    def handle(self, *args, **options):
        functions = [
            'create_resources',
            'create_protein_segments',
            'create_residue_numbering_schemes',
            'create_anomalies',
        ]

        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def create_protein_sites(self):
        self.logger.info('CREATING PROTEIN SITES')
        Site.objects.get_or_create(slug='sodium_pocket', name='Sodium ion pocket')
        self.logger.info('COMPLETED PROTEIN SITES')

    def create_resources(self):
        self.logger.info('CREATING RESOURCES')
        self.logger.info('Parsing file ' + self.resource_source_file)

        with open(self.resource_source_file, "r", encoding='UTF-8') as resource_source_file:
            for row in resource_source_file:
                split_row = shlex.split(row)

                # create resource
                try:
                    defaults = {
                        'name': split_row[1],
                        'url': split_row[2]
                    }

                    wr, created = WebResource.objects.get_or_create(slug=split_row[0], defaults=defaults)

                    if created:
                        self.logger.info('Created resource ' + wr.slug)

                except:
                    self.logger.error('Failed creating resource ' + split_row[0])
                    continue

        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('COMPLETED CREATING RESOURCES')

    def create_protein_segments(self):
        self.logger.info('CREATING PROTEIN SEGMENTS')
        self.logger.info('Parsing file ' + self.segment_source_file)

        with open(self.segment_source_file, "r", encoding='UTF-8') as segment_file:
            for row in segment_file:
                split_row = shlex.split(row)

                if int(split_row[2]):
                    fully_aligned = True
                else:
                    fully_aligned = False

                if len(split_row)==6:
                    domain = split_row[5]
                else:
                    domain = None

                # create segment
                try:
                    defaults={
                        'category': split_row[1],
                        'fully_aligned': fully_aligned,
                        'name': split_row[3]
                    }

                    s, created = ProteinSegment.objects.get_or_create(slug=split_row[0], proteinfamily=split_row[4], domain=domain, defaults=defaults)
                    s.save()

                    if created:
                        self.logger.info('Created protein segment ' + s.name)
                except Exception as msg:
                    # print('Failed creating protein segment', split_row, msg)
                    self.logger.error('Failed creating protein segment', split_row)
                    # continue
        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('COMPLETED CREATING PROTEIN SEGMENTS')

    def create_residue_numbering_schemes(self):
        self.logger.info('CREATING RESIDUE NUMBERING SCHEMES')
        self.logger.info('Parsing file ' + self.residue_number_scheme_source_file)

        with open(self.residue_number_scheme_source_file, "r", encoding='UTF-8') as residue_number_scheme_source_file:
            for row in residue_number_scheme_source_file:
                split_row = shlex.split(row)

                # create scheme
                try:
                    defaults={
                        'short_name': split_row[1],
                        'name': split_row[2]
                    }
                    if len(split_row) == 4:
                        try:
                            prns = split_row[3]
                            parent = ResidueNumberingScheme.objects.get(slug=prns)
                        except ResidueNumberingScheme.DoesNotExist:
                            raise Exception('Parent scheme {} does not exist, aborting!'.format(prns))
                        defaults['parent'] = parent

                    s, created = ResidueNumberingScheme.objects.get_or_create(slug=split_row[0], defaults=defaults)

                    if created:
                        self.logger.info('Created residue numbering scheme ' + s.short_name)
                except:
                    self.logger.error('Failed creating residue numbering scheme {}'.format(split_row[0]))
                    continue

        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('COMPLETED CREATING RESIDUE NUMBERING SCHEMES')

    def create_anomalies(self):
        self.logger.info('CREATING PROTEIN ANOMALIES')

        filenames = os.listdir(self.anomaly_source_dir)
        for source_file in filenames:
            source_file_path = os.sep.join([self.anomaly_source_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Parsing file {}'.format(source_file_path))
                # read the yaml file
                with open(source_file_path, 'r') as f:
                    ano = yaml.load(f, Loader=yaml.FullLoader)

                    # anomaly type
                    if 'anomaly_type' in ano and ano['anomaly_type']:
                        at, created = ProteinAnomalyType.objects.get_or_create(slug=ano['anomaly_type'],
                            defaults={'name': ano['anomaly_type'].title()})
                        if created:
                            self.logger.info('Created protein anomaly type {}'.format(at.slug))
                    else:
                        self.logger.error('Anomaly type not specified in file {}, skipping!'.format(source_file))
                        continue

                    # protein segment
                    if 'protein_segment' in ano and ano['protein_segment']:
                        try:
                            ps = ProteinSegment.objects.get(slug=ano['protein_segment'])
                        except ProteinSegment.DoesNotExist:
                            self.logger.error('Protein segment {} not found, skipping!'.format(
                                ano['protein_segment']))
                            continue

                    # generic number
                    if 'generic_number' in ano and ano['generic_number']:
                        rns = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
                        gn, created = ResidueGenericNumber.objects.get_or_create(label=ano['generic_number'],
                            scheme=rns, defaults={'protein_segment': ps})
                        if created:
                            self.logger.info('Created generic number {}'.format(gn.label))
                    else:
                        self.logger.error('Generic number not specified in file {}, skipping!'.format(source_file))
                        continue

                    # anomaly
                    pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=at, generic_number=gn)
                    if created:
                        self.logger.info('Created {} {}'.format(at.slug, gn.label))

                    # rule sets
                    if 'rule_sets' not in ano:
                        self.logger.error('No rule sets specified in file {}, skipping!'.format(source_file))
                        continue
                    for i, rule_set in enumerate(ano['rule_sets']):
                        # exclusive rule set?
                        exclusive = False
                        if 'exclusive' in rule_set and rule_set['exclusive']:
                            exclusive = True

                        # rules in this rule set
                        if 'rules' in rule_set and rule_set['rules']:
                            pars = ProteinAnomalyRuleSet.objects.create(protein_anomaly=pa, exclusive=exclusive)
                            self.logger.info('Created protein anomaly rule set')
                            for rule in rule_set['rules']:
                                exp_keys = ['generic_number', 'amino_acid']

                                # are all expected values specified for this rule
                                if all(x in rule for x in exp_keys) and all(rule[x] for x in exp_keys):
                                    # is this a negative rule (!= instead of ==)
                                    negative = False
                                    if rule['negative']:
                                        negative = True
                                    pargn, created = ResidueGenericNumber.objects.get_or_create(
                                        label=rule['generic_number'], protein_segment=ps, scheme=rns)
                                    if created:
                                        self.logger.info('Created generic_number {}'.format(pargn.label))
                                    par = ProteinAnomalyRule.objects.create(generic_number=pargn,
                                        rule_set=pars, amino_acid=rule['amino_acid'], negative=negative)
                                    self.logger.info('Created rule {} = {}, negative = {}'.format(
                                        par.generic_number.label, par.amino_acid, negative))
                                else:
                                    self.logger.error('Missing values for rule {:d} in file {}'.format((i+1),
                                        source_file))
        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('COMPLETED CREATING PROTEIN ANOMALIES')
