# Generated by Django 3.1.7 on 2022-05-24 10:50

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0016_experimentaldata'),
    ]

    operations = [
        migrations.AlterField(
            model_name='experimentaldata',
            name='activity_ranges',
            field=models.CharField(max_length=40, null=True),
        ),
        migrations.AlterField(
            model_name='experimentaldata',
            name='document_chembl_id',
            field=models.CharField(max_length=100, null=True),
        ),
        migrations.AlterField(
            model_name='experimentaldata',
            name='source',
            field=models.CharField(max_length=50),
        ),
        migrations.AlterField(
            model_name='experimentaldata',
            name='value_type',
            field=models.CharField(max_length=50),
        ),
    ]
