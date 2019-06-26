# Generated by Django 2.0.13 on 2019-06-21 12:05

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0005_proteingproteinpair_g_protein_subunit'),
        ('common', '0001_initial'),
        ('ligand', '0028_auto_20190621_0904'),
    ]

    operations = [
        migrations.CreateModel(
            name='AnalyzedAssay',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('order_no', models.IntegerField(max_length=2)),
                ('signalling_protein', models.CharField(max_length=40)),
                ('cell_line', models.CharField(max_length=20, null=True)),
                ('assay_type', models.CharField(max_length=50, null=True)),
                ('assay_measure', models.CharField(max_length=51, null=True)),
                ('assay_time_resolved', models.CharField(max_length=52, null=True)),
                ('ligand_function', models.CharField(max_length=53, null=True)),
                ('quantitive_measure_type', models.CharField(max_length=20, null=True)),
                ('quantitive_activity', models.FloatField(max_length=10, null=True)),
                ('quantitive_activity_sign', models.CharField(max_length=3, null=True)),
                ('quantitive_unit', models.CharField(max_length=10, null=True)),
                ('qualitative_activity', models.CharField(max_length=30, null=True)),
                ('quantitive_efficacy', models.FloatField(max_length=20, null=True)),
                ('efficacy_measure_type', models.CharField(max_length=30, null=True)),
                ('efficacy_sign', models.CharField(max_length=3, null=True)),
                ('efficacy_unit', models.CharField(max_length=20, null=True)),
                ('bias', models.CharField(max_length=5, null=True)),
                ('t_coefficient', models.FloatField(max_length=10, null=True)),
                ('t_value', models.FloatField(max_length=10, null=True)),
                ('log_bias_factor', models.FloatField(max_length=10, null=True)),
                ('emax_ligand_reference', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='ExperimentAssay.bias_ligand_reference+', to='ligand.Ligand')),
            ],
        ),
        migrations.CreateModel(
            name='AnalyzedExperiment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('submission_name', models.CharField(max_length=100)),
                ('ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.Ligand')),
                ('publication', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='common.Publication')),
                ('receptor', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='protein.Protein')),
            ],
        ),
        migrations.RemoveField(
            model_name='referenceligand',
            name='emax_ligand_reference',
        ),
        migrations.RemoveField(
            model_name='referenceligand',
            name='ligand',
        ),
        migrations.AlterField(
            model_name='experimentassay',
            name='biased_experiment',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='experiment_data', to='ligand.BiasedExperiment'),
        ),
        migrations.DeleteModel(
            name='ReferenceLigand',
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='experiment',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='analyzed_data', to='ligand.AnalyzedExperiment'),
        ),
    ]
