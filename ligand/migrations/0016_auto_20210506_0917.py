# Generated by Django 3.1.6 on 2021-05-06 07:17

from django.db import migrations, models



class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0015_experimentassay_effector_family'),
    ]

    operations = [
        migrations.AddField(
            model_name='analyzedassay',
            name='effector_family',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='measured_biological_process',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='measured_effector',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='signal_detection_tecnique',
            field=models.TextField(null=True),
        ),
    ]
