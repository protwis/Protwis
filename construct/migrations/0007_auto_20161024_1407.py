# -*- coding: utf-8 -*-
# Generated by Django 1.10.1 on 2016-10-24 12:07
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('construct', '0006_construct_snakecache'),
    ]

    operations = [
        migrations.AlterField(
            model_name='chemicalconc',
            name='concentration',
            field=models.CharField(max_length=200),
        ),
    ]
