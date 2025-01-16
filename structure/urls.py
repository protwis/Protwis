from django.conf.urls import url
from django.urls import path
from structure.views import *
from structure import views
from django.conf import settings
from django.urls import path
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

urlpatterns = [
    url(r'^$', cache_page(60*60*24)(StructureBrowser.as_view()), name='structure_browser'),
    # url(r'^$', StructureBrowser.as_view(), name='structure_browser'),
    url(r'^g_protein_structure_browser$', cache_page(60*60*24)(EffectorStructureBrowser.as_view(effector='gprot')), name='g_protein_structure_browser'),
    url(r'^arrestin_structure_browser$', cache_page(60*60*24)(EffectorStructureBrowser.as_view(effector='arrestin')), name='arrestin_structure_browser'),
    path('structure_similarity_search', StructureBlastView.as_view(), name='structure-similarity-search'),
    url(r'^browser$', RedirectBrowser, name='redirect_browser'),
    url(r'^selection_convert$', ConvertStructuresToProteins, name='convert'),
    url(r'^selection_convert_model$', ConvertStructureModelsToProteins, name='convert_mod'),
    url(r'^selection_convert_signprot_model$', ConvertStructureComplexSignprotToProteins, name='convert_signprot'),
    url(r'^hommod_download$', HommodDownload, name='hommod_download'),
    url(r'^complexmod_download$', ComplexmodDownload, name='complexmod_download'),
    url(r'^template_browser', TemplateBrowser.as_view(), name='structure_browser'),
    url(r'^template_selection', TemplateTargetSelection.as_view(), name='structure_browser'),
    url(r'^template_segment_selection', TemplateSegmentSelection.as_view(), name='structure_browser'),
    url(r'^gprot_statistics$', cache_page(60*60*24)(StructureStatistics.as_view(origin='gprot')), name='structure_statistics'),
    url(r'^arrestin_statistics$', cache_page(60*60*24)(StructureStatistics.as_view(origin='arrestin')), name='structure_statistics'),
    url(r'^statistics$', cache_page(60*60*24)(StructureStatistics.as_view()), name='structure_statistics'),
    url(r'^homology_models$', cache_page(60*60*24)(ServeHomologyModels.as_view()), name='homology_models'),
    path('ligand_complex_models', LigandComplexModels.as_view(), name='ligand_complex_models'),
    url(r'^ligand_complex_models/(?P<header>[^/]+)$', LigandComplexDetails, name='ligand_complex_details'),
    url(r'^ligand_complex_models/view/(?P<modelname>[^/]+)$', ServeComplexModDiagram, name='complexmod_serve_view'),
    path('lig_complexmod_download', LigComplexmodDownload, name='lig_complexmod_download'),
    url(r'^ligand_complex_models/(?P<modelname>[^/]+)/download_lig_complex_pdb$', SingleLigComplexModelDownload, name='single_complex_model_download'),
    # url(r'^homology_models$', ServeHomologyModels.as_view(), name='homology_models'),
    url(r'^complex_models$', cache_page(60*60*24)(ServeComplexModels.as_view()), name='complex_models'),
    url(r'^arrestin_models$', cache_page(60*60*24)(ServeComplexModels.as_view(signalling_protein='af-arrestin')), name='arrestin_complex_models'),

    # url(r'^complex_models$', ServeComplexModels.as_view(), name='complex_models'),
    url(r'^model_statistics$', cache_page(60*60*24)(ServeModelStatistics.as_view()), name='model_statistics'),
    url(r'^pdb_segment_selection', PDBSegmentSelection.as_view(), name='pdb_download'),
    url(r'^pdb_download$', PDBClean.as_view(), name='pdb_download'),
    url(r'^pdb_download_custom$', PDBClean.as_view(), name='pdb_download_custom'),
    url(r'^pdb_download/(?P<substructure>\w+)$', PDBDownload.as_view(), name='pdb_download'),
    url(r'^generic_numbering_index', GenericNumberingIndex.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results$', GenericNumberingResults.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results/(?P<substructure>\w+)$', GenericNumberingDownload.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_selection', GenericNumberingSelection.as_view(), name='generic_numbering'),
    url(r'^superposition_workflow_index$', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_gprot_index$', SuperpositionWorkflowIndex.as_view(website='gprot'), name='superposition_workflow_gprot'),
    url(r'^superposition_workflow_arrestin_index$', SuperpositionWorkflowIndex.as_view(website='arrestin'), name='superposition_workflow_gprot'),
    url(r'^segmentselectiongprot$', SegmentSelectionGprotein.as_view(), name='superposition_workflow_gprot'),
    url(r'^superposition_workflow_index/(?P<clear>\w{4})$', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_gprot_index/(?P<clear>\w{4})$', SuperpositionWorkflowIndex.as_view(website='gprot'), name='superposition_workflow_gprot'),
    url(r'^superposition_workflow_arrestin_index/(?P<clear>\w{4})$', SuperpositionWorkflowIndex.as_view(website='arrestin'), name='superposition_workflow_gprot'),
    url(r'^superposition_workflow_selection', SuperpositionWorkflowSelection.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_results$', SuperpositionWorkflowResults.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_results_gprot$', SuperpositionWorkflowResults.as_view(website='gprot'), name='superposition_workflow'),
    url(r'^superposition_workflow_results_arrestin$', SuperpositionWorkflowResults.as_view(website='arrestin'), name='superposition_workflow'),
    url(r'^superposition_workflow_results/(?P<substructure>\w+)$', SuperpositionWorkflowDownload.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_results_gprot/(?P<substructure>\w+)$', SuperpositionWorkflowDownload.as_view(website='gprot'), name='superposition_workflow'),
    url(r'^output/(?P<outfile>\w+.\w{3})/(?P<replacement_tag>\w+)$', ServePdbOutfile, name='structural_tools_result'),
    url(r'^zipoutput/(?P<outfile>\w+.\w{3})/', ServeZipOutfile, name='structural_tools_result'),
    url(r'^showtrees', RenderTrees, name='render'),
    url(r'^(?P<pdbname>\w+)$', cache_page(60*60*24*7)(StructureDetails), name='structure_details'),
    url(r'^pdb/(?P<pdbname>\w+)$', cache_page(60*60*24*7)(ServePdbDiagram), name='structure_serve_pdb'),
    url(r'^pdb_upright/(?P<pdbname>\w+)$', ServeUprightPdbDiagram, name='structure_serve_upright_pdb'),
    url(r'^pdb_clean/(?P<pdbname>\w+)_(?P<ligname>\w+)$', ServeCleanPdbDiagram, name='structure_clean_pdb'),
    url(r'^refined/(?P<modelname>\w+)_(?P<fullness>\w+)/download_pdb$', SingleModelDownload, name='single_model_download'),
    url(r'^(?P<pdbcode>\w+)/download_pdb$', SingleStructureDownload, name='single_structure_download'),
    url(r'^refined/(?P<pdbname>\w+)$', cache_page(60*60*24*7)(RefinedModelDetails), name="refined_model_details"),
    url(r'^refined/(?P<modelname>\w+)/download_complex_pdb$', SingleComplexModelDownload, name='single_complex_model_download'),
    url(r'^pdb/(?P<pdbname>\w+)/ligand/(?P<ligand>.+)$', ServePdbLigandDiagram, name='structure_serve_pdb_ligand'),
    url(r'^complex_models/(?P<header>\w+)$',cache_page(60*60*24*7)(ComplexModelDetails), name='complex_model_details'),
    url(r'^complex_models/view/(?P<modelname>\w+)$', ServeComplexModDiagram, name='complexmod_serve_view'),
]
