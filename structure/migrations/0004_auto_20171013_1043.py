# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-10-13 08:43
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('protein', '0002_auto_20170908_0758'),
        ('structure', '0003_structurerefinedseqsim_structurerefinedstatsrotamer'),
    ]

    operations = [
        migrations.CreateModel(
            name='IdentifiedSites',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.ProteinConformation')),
                ('residues', models.ManyToManyField(to='residue.Residue')),
            ],
        ),
        migrations.CreateModel(
            name='Site',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('slug', models.CharField(max_length=20)),
                ('name', models.CharField(max_length=30)),
            ],
        ),
        migrations.AddField(
            model_name='identifiedsites',
            name='site',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='structure.Site'),
        ),
    ]
