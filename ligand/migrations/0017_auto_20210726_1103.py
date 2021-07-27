# Generated by Django 3.1.6 on 2021-07-26 09:03

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0003_citation_page_name'),
        ('ligand', '0016_auto_20210722_0845'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='gtp_endogenous_ligand',
            name='publication',
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='web_links',
            field=models.ManyToManyField(to='common.WebLink'),
        ),
    ]
