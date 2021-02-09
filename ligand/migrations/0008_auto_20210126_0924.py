# Generated by Django 3.0.8 on 2021-01-26 08:24

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0032_structure_author_state'),
        ('protein', '0009_auto_20200511_1818'),
        ('ligand', '0007_ligand_pdbe'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='analyzedexperiment',
            name='chembl',
        ),
        migrations.RemoveField(
            model_name='analyzedexperiment',
            name='mutation',
        ),
        migrations.RemoveField(
            model_name='analyzedexperiment',
            name='residue',
        ),
        migrations.RemoveField(
            model_name='assayexperiment',
            name='chembl',
        ),
        migrations.RemoveField(
            model_name='assayexperiment',
            name='smiles',
        ),
        migrations.RemoveField(
            model_name='biasedexperiment',
            name='chembl',
        ),
        migrations.RemoveField(
            model_name='biasedexperiment',
            name='mutation',
        ),
        migrations.RemoveField(
            model_name='biasedexperiment',
            name='residue',
        ),
        migrations.AddField(
            model_name='analyzedexperiment',
            name='auxiliary_protein',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='bias_auxiliary_protein', to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='biasedexperiment',
            name='auxiliary_protein',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='auxiliary_protein', to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='biasedpathways',
            name='lignad_pubchem',
            field=models.CharField(max_length=40, null=True),
        ),
        migrations.AlterField(
            model_name='assayexperiment',
            name='assay_description',
            field=models.TextField(max_length=1500),
        ),
        migrations.AlterField(
            model_name='assayexperiment',
            name='published_value',
            field=models.CharField(max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='assayexperiment',
            name='standard_relation',
            field=models.CharField(max_length=10, null=True),
        ),
        migrations.AlterField(
            model_name='assayexperiment',
            name='standard_type',
            field=models.CharField(max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='assayexperiment',
            name='standard_units',
            field=models.CharField(max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='assayexperiment',
            name='standard_value',
            field=models.DecimalField(decimal_places=1, max_digits=20, null=True),
        ),
        migrations.AlterField(
            model_name='ligandproperities',
            name='sequence',
            field=models.CharField(max_length=1500, null=True),
        ),
        migrations.CreateModel(
            name='LigandReceptorStatistics',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type', models.CharField(max_length=3, null=True)),
                ('value', models.DecimalField(decimal_places=3, max_digits=9, null=True)),
                ('ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.Ligand')),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='target_protein', to='protein.Protein')),
                ('reference_protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='reference_protein', to='protein.Protein')),
            ],
        ),
        migrations.CreateModel(
            name='LigandPeptideStructure',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('chain', models.CharField(max_length=20)),
                ('ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.Ligand')),
                ('structure', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='structure.Structure')),
            ],
            options={
                'db_table': 'ligand_peptide_structure',
            },
        ),
        migrations.CreateModel(
            name='AssayExperimentSource',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('database', models.CharField(max_length=20, null=True)),
                ('database_id', models.CharField(max_length=30, null=True)),
                ('assay', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.AssayExperiment')),
            ],
        ),
    ]