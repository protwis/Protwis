# Generated by Django 3.0.3 on 2020-11-10 07:56

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0032_structure_author_state'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='TM_all',
        ),
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='date',
        ),
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='overall_all',
        ),
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='overall_backbone',
        ),
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='pdb',
        ),
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='service',
        ),
        migrations.RemoveField(
            model_name='structuremodelrmsd',
            name='version',
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='ECL1',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='ECL2',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='H8',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='ICL1',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='ICL2',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='binding_pocket',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='main_template',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, related_name='main_template', to='structure.Structure'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='structuremodelrmsd',
            name='target_structure',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, related_name='target_structure', to='structure.Structure'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='structuremodelrmsd',
            name='TM_backbone',
            field=models.DecimalField(decimal_places=1, max_digits=2, null=True),
        ),
        migrations.AlterField(
            model_name='structuremodelrmsd',
            name='homology_model',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='structure.StructureModel'),
        ),
    ]