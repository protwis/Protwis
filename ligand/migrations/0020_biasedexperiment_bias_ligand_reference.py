# Generated by Django 2.0.13 on 2019-05-31 07:07

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0019_remove_biasedexperiment_potency_ratio'),
    ]

    operations = [
        migrations.AddField(
            model_name='biasedexperiment',
            name='bias_ligand_reference',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='ExperimentAssay.bias_ligand_reference+', to='ligand.Ligand'),
        ),
    ]
