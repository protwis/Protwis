# Generated by Django 3.1.7 on 2022-09-15 07:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0005_auto_20210725_1110'),
    ]

    operations = [
        migrations.AddField(
            model_name='releasestatistics',
            name='database',
            field=models.TextField(null=True),
        ),
    ]
