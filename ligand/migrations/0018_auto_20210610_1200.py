# Generated by Django 3.1.6 on 2021-06-10 10:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0017_analyzedassay_log_bias_factor_e'),
    ]

    operations = [
        migrations.AlterField(
            model_name='biasedexperiment',
            name='auxiliary_protein',
            field=models.TextField(null=True),
        ),
    ]
