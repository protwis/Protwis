# Generated by Django 3.1.7 on 2021-12-06 07:53

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0019_auto_20211203_1404'),
    ]

    operations = [
        migrations.AddField(
            model_name='biaseddata',
            name='measured_process',
            field=models.CharField(max_length=60, null=True),
        ),
    ]
