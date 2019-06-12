# Generated by Django 2.0.13 on 2019-05-31 07:24

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0020_biasedexperiment_bias_ligand_reference'),
    ]

    operations = [
        migrations.AddField(
            model_name='experimentassay',
            name='bias_reference',
            field=models.CharField(max_length=4, null=True),
        ),
        migrations.AlterField(
            model_name='biasedexperiment',
            name='bias_assay',
            field=models.CharField(max_length=5),
        ),
        migrations.AlterField(
            model_name='biasedexperiment',
            name='bias_value',
            field=models.CharField(max_length=50, null=True),
        ),
    ]
