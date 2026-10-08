"""Focused checks for leakage, mapping, coordinate precedence and null arithmetic."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import savemat
from core.data_loader import load_session
from .common import bh_q, load_object_coords, session_phases
from .decoding import bin_phase, prepare_folds, fold_counts
from .object_cells import analyse_phase
from .run_analysis import parser


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
