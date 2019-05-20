# Generated by Django 2.0.13 on 2019-04-25 11:42

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0012_auto_20190424_1111'),
    ]

    operations = [
        migrations.AlterField(
            model_name='biasedexperiment',
            name='bias_pathway_relationship',
            field=models.CharField(max_length=90),
        ),
        migrations.AlterField(
            model_name='experimentassay',
            name='efficacy_unit',
            field=models.CharField(max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='experimentassay',
            name='quantitive_activity_sign',
            field=models.CharField(max_length=3, null=True),
        ),
        migrations.AlterField(
            model_name='experimentassay',
            name='quantitive_efficacy',
            field=models.FloatField(max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='experimentassay',
            name='quantitive_measure_type',
            field=models.CharField(max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name='experimentassay',
            name='quantitive_unit',
            field=models.CharField(max_length=20, null=True),
        ),
    ]
