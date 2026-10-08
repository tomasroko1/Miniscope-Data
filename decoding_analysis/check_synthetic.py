"""Assert independent synthetic accuracy, recall, false positives and displacement."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def check(data_dir, results_dir, min_recall=.7, max_fpr=.2):
    data = Path(data_dir); root = Path(results_dir)
    truth = pd.read_csv(data/'ground_truth.csv')
    dec = pd.read_csv(root/'decoding'/'decoding_summary.csv')
    objs = pd.read_csv(root/'object_cells'/'object_cells_summary.csv')
    disp = pd.read_csv(root/'object_cells'/'object_displacement_summary.csv')
    failures = pd.read_csv(root/'analysis_failures.csv')
    assert failures.empty, failures.to_string()
    expected_sessions = truth[['animal','session']].drop_duplicates()
    assert len(dec) == len(expected_sessions)*4, 'faltan fases/sesiones de decoding'
    ok = (dec.median_error_cm < .85*dec.shuffle_p5_cm) & (dec.median_error_cm < .85*dec.chance_median_error_cm)
    print(dec.groupby('phase')[['median_error_cm','shuffle_p5_cm','chance_median_error_cm']].median().round(2).to_string())
    assert ok.all(), 'Decoding no supera claramente nulo/azar:\n'+dec[~ok].to_string()
    for _, group in objs.groupby(['animal', 'session_file']):
        assert set(group.phase) == {'OF1', 'SAMPLE', 'TEST', 'OF2'}, 'faltan controles/fases de objetos'
    assert not objs.coords_are_placeholder.any(), 'coordenadas sintéticas marcadas placeholder'
    sample = objs[objs.phase == 'SAMPLE']
    merged = truth.merge(sample, left_on=['animal','session','global_id'],
                         right_on=['animal','session_file','global_id'], how='inner')
    planted = merged[merged.cell_type == 'object']; others = merged[merged.cell_type != 'object']
    assert len(planted) == (truth.cell_type == 'object').sum(), 'faltan células plantadas'
    recall = planted.is_object_cell.mean(); fpr = others.is_object_cell.mean()
    print('Recall SAMPLE: %.3f (%d/%d); falsos positivos: %.3f (%d/%d)' %
          (recall,planted.is_object_cell.sum(),len(planted),fpr,others.is_object_cell.sum(),len(others)))
    assert recall >= min_recall, 'recall insuficiente'
    assert fpr <= max_fpr, 'falsos positivos excesivos'
    aligned = planted[planted.is_object_cell].merge(disp, on=['animal','session_file','global_id'], suffixes=('_sample','_disp'))
    assert len(aligned) == planted.is_object_cell.sum(), 'faltan resultados de desplazamiento'
    assert (aligned.behaviour_in_TEST == aligned.displacement).all(), 'clasificación desplazamiento incorrecta'
    moved = aligned[aligned.displacement == 'moved']; stayed = aligned[aligned.displacement == 'stayed']
    assert len(moved) and len(stayed), 'faltan ejemplos moved/stayed'
    assert (moved.follow_index > 0).mean() >= .8, 'células moved no siguen objeto'
    print('Desplazamiento correcto: moved=%d, stayed=%d; follow_index mediano moved=%.3f; significativas TEST=%.3f' %
          (len(moved),len(stayed),moved.follow_index.median(),aligned.significant_in_TEST.mean()))
    print('PASS: todas las comprobaciones sintéticas.')


def main():
    p = argparse.ArgumentParser(); p.add_argument('--data-dir',required=True); p.add_argument('--results-dir',required=True)
    p.add_argument('--min-recall',type=float,default=.7); p.add_argument('--max-fpr',type=float,default=.2)
    a=p.parse_args(); check(a.data_dir,a.results_dir,a.min_recall,a.max_fpr)


if __name__ == '__main__': main()
