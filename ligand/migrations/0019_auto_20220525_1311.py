# Generated by Django 3.1.7 on 2022-05-25 11:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0018_auto_20220524_1654'),
    ]

    operations = [
        migrations.RenameField(
            model_name='experimentaldata',
            old_name='activity_ranges',
            new_name='p_activity_ranges',
        ),
        migrations.RenameField(
            model_name='experimentaldata',
            old_name='activity_value',
            new_name='p_activity_value',
        ),
        migrations.AddField(
            model_name='experimentaldata',
            name='standard_activity_value',
            field=models.CharField(max_length=10, null=True),
        ),
    ]
