# Generated by Django 3.1.6 on 2021-05-12 18:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0017_auto_20210511_0135'),
    ]

    operations = [
        migrations.AlterField(
            model_name='analyzedassay',
            name='assay_description',
            field=models.CharField(max_length=901, null=True),
        ),
    ]
