# Generated by Django 3.1.6 on 2021-05-12 21:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0016_auto_20210506_0917'),
    ]

    operations = [
        migrations.AddField(
            model_name='analyzedassay',
            name='log_bias_factor_a',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='log_bias_factor_b',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='log_bias_factor_c',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='log_bias_factor_d',
            field=models.CharField(max_length=60, null=True),
        ),
    ]
