# Generated by Django 3.1.7 on 2022-10-12 08:35

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0021_auto_20221005_1239'),
    ]

    operations = [
        migrations.AddField(
            model_name='assayexperiment',
            name='reference_ligand',
            field=models.CharField(max_length=300, null=True),
        ),
    ]
