# Generated by Django 2.0.13 on 2019-06-25 20:16

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0035_auto_20190625_1406'),
    ]

    operations = [
        migrations.AlterField(
            model_name='analyzedassay',
            name='potency',
            field=models.CharField(max_length=10, null=True),
        ),
    ]
