# Generated by Django 3.0.8 on 2020-10-28 08:08

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0010_assayexperiment_dublicate'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='assayexperiment_dublicate',
            name='chembl',
        ),
        migrations.RemoveField(
            model_name='assayexperiment_dublicate',
            name='published_relation',
        ),
        migrations.RemoveField(
            model_name='assayexperiment_dublicate',
            name='published_type',
        ),
        migrations.RemoveField(
            model_name='assayexperiment_dublicate',
            name='published_units',
        ),
        migrations.RemoveField(
            model_name='assayexperiment_dublicate',
            name='published_value',
        ),
        migrations.RemoveField(
            model_name='assayexperiment_dublicate',
            name='smiles',
        ),
    ]