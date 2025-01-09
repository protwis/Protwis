from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db.models import F, Q
from common.tools import test_model_updates
from common.models import ReleaseNotes, ReleaseStatistics, ReleaseStatisticsType
from drugs.models import Drugs
from residue.models import ResidueGenericNumber
from interaction.models import ResidueFragmentInteraction
from ligand.models import Ligand, AssayExperiment, BiasedData, BiasedPathwaysAssay, Endogenous_GTP, BalancedLigands
from mutation.models import MutationExperiment
from mutational_landscape.models import NaturalMutations
from protein.models import Protein, ProteinCouplings
from structure.models import Structure, StructureModel, StructureComplexModel
from signprot.models import SignprotComplex, SignprotStructure
from contactnetwork.models import InteractingResiduePair

import logging
import shlex
import os
import yaml
import django.apps

class Command(BaseCommand):
    help = 'Builds release notes and stats'

    logger = logging.getLogger(__name__)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    release_notes_data_dir = os.sep.join([settings.DATA_DIR, 'release_notes'])

    def handle(self, *args, **options):
        try:
            self.create_release_notes()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
        test_model_updates(self.all_models, self.tracker, check=True)

    def create_release_notes(self):
        self.logger.info('CREATING RELEASE NOTES')
        ReleaseNotes.objects.all().delete()

        # what files should be parsed?
        filenames = os.listdir(self.release_notes_data_dir)

        for source_file in filenames:
            # self.logger.info('Parsing file ' + source_file)
            source_file_path = os.sep.join([self.release_notes_data_dir, source_file])

            if source_file.endswith('.html'):
                file_date = source_file_path[:-5].split("/")[-1]
                release_notes, created = ReleaseNotes.objects.get_or_create(date=file_date)
                # if created:
                    # self.logger.info('Created release notes for {}'.format(file_date))

                with open(source_file_path[:-4]+'html', 'r') as h:
                    release_notes.html = h.read()
                    release_notes.save()

                    if created:
                        self.logger.info('Updated html for release notes {}'.format(file_date))

        # generate statistics
        latest_release_notes = ReleaseNotes.objects.all()[0]

        # G protein structures (both complexes and individual)
        complexstructs = list(SignprotComplex.objects.filter(protein__family__slug__startswith='100').values_list("structure__pdb_code__index", flat = True))
        ncstructs = list(SignprotStructure.objects.filter(protein__family__slug__startswith='100').values_list("pdb_code__index", flat = True))
        gprotein_structs = set(complexstructs + ncstructs)
        signcomp = SignprotComplex.objects.all().exclude(protein__family__slug__startswith="200")
        interface_interactions_count = (
            InteractingResiduePair.objects.filter(referenced_structure__in=signcomp.values_list("structure", flat=True))
            .exclude(res1__protein_conformation_id=F("res2__protein_conformation_id"))
            .count()
        )
        signcomp_arrestin = SignprotComplex.objects.filter(protein__family__slug__startswith="200")
        interface_interactions_arrestin_count = (
            InteractingResiduePair.objects.filter(referenced_structure__in=signcomp_arrestin.values_list("structure", flat=True))
            .exclude(res1__protein_conformation_id=F("res2__protein_conformation_id"))
            .count()
        )
            #['Proteins', Protein.objects.filter(sequence_type__slug='wt').count(), 'GPCRdb'],
            #['Species', Species.objects.all().count(), 'GPCRdb'],
            #['Exp. Gprotein structures', Structure.objects.filter(protein_conformation__protein__family__slug__startswith="100").count()],
            # ['Exp. Gprotein structures', len(gprotein_structs), 'GPCRdb'],
            #['Exp. Arrestin structures', Structure.objects.filter(protein_conformation__protein__family__slug__startswith="200").count()],
            # ['GPCR-G protein structure models', StructureComplexModel.objects.filter(receptor_protein__accession__isnull=False).count(), 'GPCRdb'],
        stats = [
            #GPCRdb Block
            ['Human proteins GPCRdb', Protein.objects.filter(sequence_type__slug='wt', family__slug__startswith='0', species__common_name="Human").count(), 'GPCRdb'],
            ['Species orthologs GPCRdb', Protein.objects.filter(sequence_type__slug='wt', family__slug__startswith='0').exclude(species__common_name="Human").count(), 'GPCRdb'],
            ['Generic residues GPCRdb', ResidueGenericNumber.objects.filter(scheme_id__in=[7,8,9,10,11]).values('label').count(), 'GPCRdb'],
            ['GPCR structures GPCRdb', Structure.objects.filter(protein_conformation__protein__family__slug__startswith="0").exclude(structure_type__slug__startswith='af-').count(), 'GPCRdb'],
            ['Refined structures GPCRdb', StructureModel.objects.filter(protein__accession__isnull=True, protein__family__slug__startswith="0").count() + Structure.objects.filter(structure_type__slug__startswith="af-signprot-refined").count(), 'GPCRdb'],
            ['GPCR structure models GPCRdb', StructureModel.objects.filter(protein__accession__isnull=False).count(), 'GPCRdb'],
            ['Physiological ligand-GPCR structure models GPCRdb', Structure.objects.filter(structure_type__slug__in=['af-rfaa-sm','af-signprot-peptide']).count(), 'GPCRdb'],
            ['Drugs GPCRdb', Drugs.objects.values('name').distinct().count(), 'GPCRdb'],
            ['Drug targets GPCRdb', Drugs.objects.values('target').distinct().count(), 'GPCRdb'],
            ['Disease indications GPCRdb', Drugs.objects.values('indication').distinct().count(), 'GPCRdb'],
            ['Ligands GPCRdb', Ligand.objects.all().count(), 'GPCRdb'],
            ['Physiological ligands GPCRdb', Endogenous_GTP.objects.values('ligand_id').distinct().count(), 'GPCRdb'],
            ['Ligand bioactivities GPCRdb', AssayExperiment.objects.all().count(), 'GPCRdb'],
            ['Ligand site mutations GPCRdb', MutationExperiment.objects.all().count(), 'GPCRdb'],
            ['Ligand interactions GPCRdb', ResidueFragmentInteraction.objects.all().exclude(structure_ligand_pair__structure__structure_type__slug__startswith='af-').count(), 'GPCRdb'],
            #GproteinDb block
            ['Human G proteins GproteinDb', Protein.objects.filter(family__parent__parent__name="Alpha", species__common_name="Human", accession__isnull=False).count(), 'GproteinDb'],
            ['Species orthologs GproteinDb', Protein.objects.filter(family__parent__parent__name="Alpha", accession__isnull=False).count(), 'GproteinDb'],
            ['G protein couplings GproteinDb', ProteinCouplings.objects.all().exclude(g_protein__slug__startswith="200").count(), 'GproteinDb'],
            ['G proteins GproteinDb', signcomp.exclude(structure__structure_type__slug__startswith='af-').count() + SignprotStructure.objects.all().exclude(protein__family__slug__startswith="200").count(), 'GproteinDb'],
            ['G protein complexes GproteinDb', SignprotComplex.objects.all().exclude(structure__structure_type__slug__startswith='af-').count(), 'GproteinDb'],
            ['Generic residues GproteinDb', ResidueGenericNumber.objects.filter(scheme_id__in=[15]).values('label').count(), 'GproteinDb'],
            ['GPCR complexes', Structure.objects.filter(structure_type__slug='af-signprot').count(), 'GproteinDb'],
            ['Refined complexes GproteinDb', Structure.objects.filter(structure_type__slug__startswith='af-signprot-refined').count(), 'GproteinDb'],
            ['G protein interface GproteinDb', interface_interactions_count, 'GproteinDb'],
            ['Interface mutations GproteinDb', 54, 'GproteinDb'],
            #ArrestinDb block
            ['Human arrestins ArrestinDb', Protein.objects.filter(family__slug__startswith="200", species__common_name="Human", accession__isnull=False).count(), 'ArrestinDb'],
            ['Species orthologs ArrestinDb', Protein.objects.filter(family__slug__startswith="200", accession__isnull=False).count(), 'ArrestinDb'],
            ['Arrestin couplings ArrestinDb', ProteinCouplings.objects.filter(g_protein__slug__startswith="200").count(), 'ArrestinDb'],
            ['Arrestins ArrestinDb', signcomp_arrestin.count() + SignprotStructure.objects.filter(protein__family__slug__startswith="200").count(), 'ArrestinDb'],
            ['GPCR complexes ArrestinDb', signcomp_arrestin.count(), 'ArrestinDb'],
            ['Generic residues ArrestinDb', ResidueGenericNumber.objects.filter(scheme_id__in=[16]).values('label').count(), 'ArrestinDb'],
            ['Interface interactions ArrestinDb', interface_interactions_arrestin_count, 'ArrestinDb'],
            ['Interface mutations ArrestinDb', 409, 'ArrestinDb'],
            #BiasedSignalingAtlas block
            ['Physiology-biased Biased Signaling Atlas', BiasedData.objects.filter(pathway_biased__isnull=False).values_list("ligand_id").distinct().count(), 'Biased Signaling Atlas'],
            ['Pathway-biased Biased Signaling Atlas', BiasedData.objects.filter(pathway_biased__isnull=False).values_list("ligand_id").distinct().count(), 'Biased Signaling Atlas'],
            ['Total datapoints Biased Signaling Atlas', BiasedData.objects.all().count(), 'Biased Signaling Atlas'],
            ['Pathway effects Biased Signaling Atlas', BiasedPathwaysAssay.objects.all().count(), 'Biased Signaling Atlas'],
            ['Ligands Biased Signaling Atlas', BiasedData.objects.filter(pathway_preferred__isnull=False).values_list("ligand_id").distinct().count(), 'Biased Signaling Atlas'],
            ['Total datapoints Biased Signaling Atlas', BiasedData.objects.filter(pathway_preferred__isnull=False).count(), 'Biased Signaling Atlas'],
            ['For pathway-bias Biased Signaling Atlas', BalancedLigands.objects.all().values_list("ligand_id").count(), 'Biased Signaling Atlas'],
            ['For physiology-bias Biased Signaling Atlas', Endogenous_GTP.objects.filter(Q(endogenous_status="Principal") | Q(potency_ranking=1)).count(), 'Biased Signaling Atlas']
        ]
        for stat in stats:
            stat_type, created = ReleaseStatisticsType.objects.get_or_create(name=stat[0])
            if created:
                self.logger.info('Created release statistics type {}'.format(stat[0]))

            rstat, created = ReleaseStatistics.objects.update_or_create(release=latest_release_notes, statistics_type=stat_type,
                defaults={'value': stat[1], 'database': stat[2]})
            if created:
                self.logger.info('Created release stat {} {} {}'.format(latest_release_notes, stat[1], stat[2]))


        self.logger.info('COMPLETED CREATING RELEASE NOTES')
