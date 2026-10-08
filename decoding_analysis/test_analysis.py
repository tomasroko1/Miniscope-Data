"""Focused checks for leakage, mapping, coordinate precedence and null arithmetic."""
import tempfile
import contextlib
import io
import json
import warnings
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.io import savemat
from core.data_loader import load_session
from core.place_fields import threshold_event_frames
from .common import activity, bh_q, load_object_coords, session_phases
from .decoding import bin_phase, prepare_folds, fold_counts, decode_phase, cross_validate
from .object_cells import analyse_phase, flag_of_controls, example_figures
from .run_analysis import parser, decoding_reports, main as run_main


class AnalysisTests(unittest.TestCase):
    def options(self):
        opts = parser().parse_args(['all', '--n-examples', '0'])
        opts.decode_min_speed = 2; opts.object_min_speed = 0
        opts.decode_min_events = 1; opts.object_min_events = 10
        return opts

    def phase(self, duration=80):
        t = np.arange(int(duration*20))/20
        rng = np.random.default_rng(2)
        return dict(phase='SAMPLE', t=t, dt=np.full(len(t),.05), x=15*np.sin(t),
                    y=10*np.cos(t), S=(rng.random((len(t),3)) < .025).astype(float),
                    speed=np.full(len(t),5.), valid=np.ones(len(t),bool),
                    ids=np.arange(3), columns=np.arange(3), aligned=True)

    def test_empty_and_incomplete_bins(self):
        phase = self.phase()
        phase['t'] = phase['t'].copy(); phase['t'][800:] += 3
        phase['dt'] = np.r_[np.diff(phase['t']), .05]
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            b = bin_phase(phase, .5)
        empty = np.asarray(b['aggregate'].sum(axis=1)).ravel() == 0
        self.assertEqual(empty.sum(), 6)
        self.assertFalse(b['valid'][empty].any())
        self.assertFalse(b['valid'][79])
        self.assertTrue(np.isfinite(b['xy']).all())
        self.assertTrue(b['valid'][86:].all())
        opts = self.options(); opts.decode_shuffles = 2
        metrics, _, _ = decode_phase(b, opts, np.random.default_rng(4))
        self.assertGreater(metrics['n_bins'], 0)

    def test_vector_events_equal_core(self):
        rng = np.random.default_rng(14)
        S = rng.normal(size=(701, 27)); source = rng.normal(size=(219, 27))
        S[rng.random(S.shape) < .07] = np.nan
        source[rng.random(source.shape) < .11] = np.nan
        S[0, 0] = np.inf; source[0, 0] = -np.inf
        source[:, 1] = 0; source[2:, 2] = np.nan
        for src in [None, source]:
            expected = np.column_stack([threshold_event_frames(S[:, i],
                       threshold_source=None if src is None else src[:, i])[0]
                       for i in range(S.shape[1])])
            np.testing.assert_array_equal(activity(S, source=src), expected)
        source[:, 2] = np.nan
        with self.assertRaisesRegex(ValueError, 'insuficientes datos'):
            activity(S, source=source)

    def test_arena_projection_and_gap_occupancy(self):
        phase = self.phase(); phase['x'][:] = 0; phase['y'][:] = 0
        phase['x'][:48] = 26; phase['x'][48:160] = 28
        phase['x'][160:208] = np.nan
        phase['t'][800:] += 3
        data = dict(phase, type='simple', v=phase['speed'])
        with contextlib.redirect_stdout(io.StringIO()) as output:
            _, phases = session_phases(data)
        ph = phases[0]
        self.assertIn('ARENA', output.getvalue())
        self.assertAlmostEqual(ph['frac_frames_projected'], .03)
        self.assertAlmostEqual(ph['frac_frames_outside_arena'], .07)
        self.assertTrue(ph['valid'][:48].all())
        self.assertTrue((np.hypot(ph['x'][:48], ph['y'][:48]) < 25).all())
        self.assertFalse(ph['valid'][48:208].any())
        self.assertFalse(ph['valid'][799]); self.assertAlmostEqual(ph['dt'][799], .05)
        self.assertEqual(data['x'][0], 26)
        _, no_tolerance = session_phases(data, arena_tolerance_cm=0)
        self.assertFalse(no_tolerance[0]['valid'][:48].any())
        data.pop('v')
        _, calculated_speed = session_phases(data)
        self.assertTrue(np.isfinite(calculated_speed[0]['speed']).all())

    def test_strict_control_flags(self):
        table = pd.DataFrame([dict(phase=p, global_id=g, is_object_cell=sig)
                              for p, g, sig in [('SAMPLE', 0, True), ('OF1', 0, False),
                              ('OF2', 0, False), ('SAMPLE', 1, True), ('OF1', 1, True),
                              ('OF2', 1, False), ('SAMPLE', 2, True), ('OF1', 2, False)]])
        flagged = flag_of_controls(table.copy(), True).query("phase == 'SAMPLE'")
        self.assertEqual(list(flagged.is_object_cell_strict), [True, False, False])
        self.assertTrue(pd.isna(flagged.sig_in_OF_control.iloc[2]))
        unavailable = flag_of_controls(table.copy(), False)
        self.assertTrue(unavailable.sig_in_OF_control.isna().all())
        self.assertFalse(unavailable.is_object_cell_strict.any())
        missing = flag_of_controls(table[table.phase != 'OF2'].copy(), True)
        self.assertFalse(missing.is_object_cell_strict.any())

    def test_shuffle_cache_and_used_cells(self):
        opts = self.options(); opts.decode_shuffles = 8
        b = bin_phase(self.phase(duration=65), opts.time_bin)
        seed = 10
        metrics, prediction, null = decode_phase(b, opts, np.random.default_rng(seed))
        folds, pixels = prepare_folds(b, opts, b['valid'])
        rng = np.random.default_rng(seed)
        cross_validate(b, folds, pixels, opts, rng=rng)
        shifts = rng.integers(60, len(b['xy'])-60+1, size=8)
        direct = [np.median(cross_validate(b, folds, pixels, opts, shift=int(s))) for s in shifts]
        np.testing.assert_array_equal(null, direct)
        self.assertEqual(metrics['n_cells_used'], prediction.groupby('fold').n_cells_fold.first().median())
        opts.activity = 'S'
        before = fold_counts(b, folds[0]['train'], opts, 61)
        opts.decode_shuffles = 0
        decode_phase(b, opts, np.random.default_rng(seed))
        np.testing.assert_array_equal(before, fold_counts(b, folds[0]['train'], opts, 61))

    def test_example_nearly_constant_null(self):
        phase = self.phase(); opts = self.options(); opts.n_examples = 1
        opts.object_shuffles = 5
        coords = pd.DataFrame([dict(object_id='A1', x_cm=0, y_cm=0, is_placeholder=0)])
        table, cache = analyse_phase(phase, coords, opts, np.random.default_rng(4))
        table = flag_of_controls(table, False)
        cache['speed_mask'] = np.ones(len(phase['t']), bool)
        cache['null'][:] = -1
        cache['null'][0] = np.nextafter(-1., 0.)
        with tempfile.TemporaryDirectory() as tmp:
            example_figures(table, {'SAMPLE': cache}, Path(tmp), opts)
            self.assertEqual(len(list(Path(tmp).glob('example_*.png'))), 1)

    def test_comparison_shuffle_z(self):
        rows = [dict(animal='M', session_file='session.mat', treatment='VEH', task='SD',
                     phase=p, median_error_cm=e, shuffle_z=z,
                     shuffle_median_of_medians_cm=20, chance_median_error_cm=20)
                for p, e, z in [('OF1', 5, -2), ('OF2', 7, -4),
                                ('SAMPLE', 4, -6), ('TEST', 6, -8)]]
        with tempfile.TemporaryDirectory() as tmp:
            decoding_reports(pd.DataFrame(rows), Path(tmp))
            delta = pd.read_csv(Path(tmp)/'decoding_objects_vs_of.csv').iloc[0]
            self.assertEqual(delta.delta_cm, -1)
            self.assertEqual(delta.delta_shuffle_z, -4)

    def test_cli_session_suffix_and_provenance(self):
        phase = self.phase()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp)/'data'/'M'; data.mkdir(parents=True)
            savemat(data/'recording.day.mat', dict(pos=dict(x=phase['x'], y=phase['y'],
                    v=phase['speed']), act=dict(S=phase['S'], t=phase['t']),
                    sess=dict(day_name='recording')))
            for i, session in enumerate(['recording.day', 'recording.day.mat']):
                out = Path(tmp)/('results%d' % i)
                args = ['run_analysis', 'decode', '--data-dir', str(data.parent),
                        '--out-dir', str(out), '--session', session, '--n-shuffles', '0',
                        '--arena-tolerance-cm', '1.5', '--n-examples', '0']
                with patch('sys.argv', args), contextlib.redirect_stdout(io.StringIO()):
                    run_main()
                config = json.loads((out/'run_config.json').read_text())
                self.assertEqual(config['session'], session)
                self.assertEqual(config['arena_tolerance_cm'], 1.5)
                self.assertEqual(config['decode_shuffles'], 0)
                self.assertEqual(config['object_shuffles'], 0)
                summary = pd.read_csv(out/'decoding'/'decoding_summary.csv')
                self.assertEqual(summary.session_file.iloc[0], 'recording.day.mat')
                self.assertEqual(summary.frac_frames_projected.iloc[0], 0)
                self.assertIn('n_cells_used', summary)

    def test_training_never_uses_test_or_gap(self):
        opts = self.options(); binned = bin_phase(self.phase(), opts.time_bin)
        folds, _ = prepare_folds(binned,opts,binned['valid'])
        for fold in folds:
            counts = fold_counts(binned,fold['train'],opts,0)
            altered = dict(binned, S=binned['S'].copy())
            not_training = ~fold['train'][np.minimum(binned['labels'], len(fold['train'])-1)]
            altered['S'][not_training] = 1e6
            after = fold_counts(altered,fold['train'],opts,0)
            np.testing.assert_array_equal(counts[fold['train']],after[fold['train']])
            self.assertFalse(fold['train'][fold['test']].any())
            distance = np.abs(binned['times'][fold['train'],None]-binned['times'][fold['test']][None,:])
            self.assertTrue((distance > opts.fold_gap_s).all())
            self.assertEqual(fold['hist'][:, fold['test']].nnz,0)

    def test_shuffle_threshold_uses_shifted_training_source(self):
        opts = self.options(); b = bin_phase(self.phase(),opts.time_bin)
        train = np.arange(len(b['xy'])) < 50; shift = 65
        before = fold_counts(b,train,opts,shift)
        changed = dict(b, S=b['S'].copy())
        source = np.roll(train,-shift)
        excluded = ~source[np.minimum(b['labels'],len(source)-1)]
        changed['S'][excluded] = 1e7
        after = fold_counts(changed,train,opts,shift)
        np.testing.assert_array_equal(before[train],after[train])

    def test_fft_null_matches_direct_circular_sums(self):
        phase = self.phase(); phase['x'][:400] = 0; phase['y'][:400] = 0
        coords = pd.DataFrame([dict(object_id='A1',x_cm=0,y_cm=0,is_placeholder=0)])
        opts = self.options(); opts.object_shuffles=8
        rng = np.random.default_rng(3)
        _, cache = analyse_phase(phase,coords,opts,rng)
        expected_shifts = np.random.default_rng(3).integers(600,len(phase['x'])-600+1,size=8)
        inside = phase['x']**2+phase['y']**2 <= opts.radius**2
        oi = phase['dt'][inside].sum(); oo=phase['dt'][~inside].sum()
        direct=[]
        for shift in expected_shifts:
            shifted = np.roll(cache['values'],-int(shift),axis=0)
            rin=shifted[inside].sum(axis=0)/oi; rout=shifted[~inside].sum(axis=0)/oo
            direct.append((rin-rout)/(rin+rout))
        np.testing.assert_allclose(cache['null'],direct,atol=1e-12)

    def test_bh_and_missing_occupancy(self):
        np.testing.assert_allclose(bh_q([.01,.04,np.nan,.03]),[.03,.04,np.nan,.04])
        phase=self.phase(); phase['x'][:]=20; phase['y'][:]=0
        coords=pd.DataFrame([dict(object_id='A1',x_cm=0,y_cm=0,is_placeholder=0)])
        opts=self.options(); opts.object_shuffles=2
        table,_=analyse_phase(phase,coords,opts,np.random.default_rng(4))
        self.assertTrue(table.p_value.isna().all())
        self.assertTrue((table.exclusion_reason == 'ocupación insuficiente').all())

    def test_coordinates_select_complete_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=[]
            for animal,day,ids in [('', '', ['A1','A2']),('', 'T1_SD_CNO',['D1']),('M','T1_SD_CNO',['M1'])]:
                for oid in ids:
                    rows.append(dict(session_type='SD',animal=animal,day_name=day,phase='SAMPLE',
                                     object_id=oid,x_cm=0,y_cm=0,match_id='',is_placeholder=0,note='test'))
            path=Path(tmp)/'coords.csv'; pd.DataFrame(rows).to_csv(path,index=False)
            self.assertEqual(list(load_object_coords(path,'M','T1_SD_CNO','SAMPLE').object_id),['M1'])
            self.assertEqual(list(load_object_coords(path,'N','T1_SD_CNO','SAMPLE').object_id),['D1'])
            self.assertEqual(list(load_object_coords(path,'N','T2_SD_VEH','SAMPLE').object_id),['A1','A2'])

    def test_simple_loader_fallback_and_ambiguous_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'M'; path.mkdir()
            savemat(path/'simple.mat',dict(pos=dict(x=np.arange(10),y=np.zeros(10)),
                    act=dict(S=np.zeros((10,2))),sess=dict(day_name='simple')))
            data=load_session('M','simple.mat',tmp)
            self.assertEqual(data['type'],'simple'); self.assertAlmostEqual(data['t'][1],.05)
            _, phases=session_phases(data,common=True)
            self.assertEqual(len(phases),1); self.assertEqual(phases[0]['S'].shape,(10,2))
            merged=dict(type='merged',subsessions=[data,data],mapping=np.zeros((2,2)),sess=None)
            _, phases=session_phases(merged,common=True)
            self.assertTrue(all(not phase['aligned'] for phase in phases))
            self.assertTrue(all(phase['S'].shape[1]==2 for phase in phases))


if __name__ == '__main__': unittest.main()
