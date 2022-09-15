# Generated by Django 3.0.3 on 2022-09-14 13:21

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0041_auto_20220801_1512'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structuremodelrmsd',
            name='main_template',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='main_template', to='structure.Structure'),
        ),
        migrations.AlterField(
            model_name='structuremodelrmsd',
            name='target_structure',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='target_structure', to='structure.Structure'),
        ),
    ]
