# Generated by Django 2.0.13 on 2019-05-06 10:38

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0014_biasedexperiment_potency_ratio'),
    ]

    operations = [
        migrations.AlterField(
            model_name='biasedexperiment',
            name='receptor',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='protein.Protein'),
        ),
    ]
