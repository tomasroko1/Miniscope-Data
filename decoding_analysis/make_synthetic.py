"""Generate MATLAB cells compatible with the real loader; no analysis shortcuts."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import savemat
from core.data_loader import load_session, mapping_columns


def geometry(task):
    separation = (15,25) if task == 'SD' else (5,15)
    sample = [('A1',-separation[0]/2,-10,''), ('A2',separation[0]/2,-10,''),
              ('A3',-separation[1]/2,10,''), ('A4',separation[1]/2,10,'')]
    test = [('T1',sample[0][1],sample[0][2],'A1'), ('T2',13,2,'A3')]
    return sample, test


def coordinate_rows(task, placeholder=0):
    sample, test = geometry(task)
    return [dict(session_type=task, day_name='', animal='', phase=phase, object_id=oid,
                 x_cm=x, y_cm=y, match_id=match, is_placeholder=placeholder,
                 note='PLACEHOLDER: sustituir por coordenadas medidas' if placeholder else 'Geometría sintética')
            for phase, coords in [('SAMPLE',sample),('TEST',test)] for oid,x,y,match in coords]


def walk(rng, duration, objects, hz=20):
    n = int(duration*hz); dt = 1/hz
    xy = np.empty((n,2)); xy[0] = [0,0]
    velocity = np.zeros(2); target = np.zeros(2); hold = 0; timer = 0
    for i in range(1,n):
        if timer <= 0:
            if objects and rng.random() < .38:
                obj = objects[rng.integers(len(objects))]
                target = np.array(obj[1:3], float) + rng.normal(0,1,2)
                hold = rng.uniform(1,3)
            else:
                angle = rng.uniform(0,2*np.pi); radius = np.sqrt(rng.random())*24
                target = radius*np.array([np.cos(angle), np.sin(angle)]); hold = rng.uniform(0,.4)
            timer = 6
            desired_speed = rng.uniform(8,22)
        delta = target-xy[i-1]; distance = np.linalg.norm(delta)
        if distance < 1.5:
            direction = rng.normal(0,.5,2)
            hold -= dt
            if hold <= 0: timer = 0
        else:
            direction = delta/max(distance,1e-9)*min(desired_speed, distance*4)
            timer -= dt
        velocity = .85*velocity+.15*direction
        xy[i] = xy[i-1]+velocity*dt
        if np.linalg.norm(xy[i]) > 24.5:
            xy[i] *= 24.5/np.linalg.norm(xy[i]); velocity *= -.5; timer = 0
    speed = np.r_[np.linalg.norm(np.diff(xy,axis=0),axis=1)/dt,0]
    return xy, speed


def generate(out, seed=0, duration=240, n_cells=80, n_animals=2):
    if duration <= 65 or n_cells < 20: raise ValueError('duration > 65 s; n-cells >= 20')
    rng = np.random.default_rng(seed); out = Path(out); out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(coordinate_rows('SD')+coordinate_rows('XsS')).to_csv(out/'object_coords_synthetic.csv', index=False)
    truth = []
    days = ['HabL','T1_SD_CNO','T2_SD_VEH','T3_XsS_CNO','T4_XsS_VEH']
    for ai in range(n_animals):
        animal = 'Synth%02d' % (ai+1); (out/animal).mkdir(exist_ok=True)
        for day in days:
            task = 'XsS' if 'XsS' in day else 'SD'
            sample, test = geometry(task)
            phases = ['OF1','OF2','OF3','OF4'] if day == 'HabL' else ['OF1','SAMPLE','TEST','OF2']
            np_place = int(.6*n_cells); np_object = min(8,n_cells-np_place) if day != 'HabL' else 0
            angles = rng.uniform(0,2*np.pi,np_place); radial = np.sqrt(rng.random(np_place))*24
            centres = np.column_stack((radial*np.cos(angles),radial*np.sin(angles)))
            sigma = rng.uniform(3,6,np_place)
            fields = {k: np.empty(4, dtype=object) for k in ['x','y','v','t','S']}
            mapping = np.full((n_cells,4),np.nan)
            for gid in range(n_cells):
                kind = 'place' if gid < np_place else ('object' if gid < np_place+np_object else 'noise')
                oid = 'A1' if (gid-np_place)%2 == 0 else 'A3'
                truth.append(dict(animal=animal, session=day+'.mat', global_id=gid, cell_type=kind,
                                  object_id=oid if kind == 'object' else '',
                                  behaviour_in_TEST=('stayed' if oid == 'A1' else 'moved') if kind == 'object' else ''))
            for pi, phase in enumerate(phases):
                objects = sample if phase == 'SAMPLE' else (test if phase == 'TEST' else [])
                xy, v = walk(rng,duration,objects); n = len(xy)
                rate = np.full((n,n_cells),.06)
                distance = ((xy[:,None,:]-centres[None,:,:])**2).sum(axis=2)
                rate[:,:np_place] += 16*np.exp(-distance/(2*sigma[None,:]**2))
                for j in range(np_object):
                    gid = np_place+j; oid = 'A1' if j%2 == 0 else 'A3'
                    if phase in ['SAMPLE','TEST']:
                        obj = next(o for o in objects if (o[0] == oid if phase == 'SAMPLE' else o[3] == oid))
                        d = np.linalg.norm(xy-np.array(obj[1:3]),axis=1)
                        rate[:,gid] = .015 + 4*(d < 4.5)
                    else: rate[:,gid] = .015
                # Mostly-zero, positive deconvolved-like events with small amplitude variability.
                S = rng.poisson(rate/20)*rng.uniform(.9,1.1,rate.shape)
                present = np.ones(n_cells,bool)
                # Keep planted cells present; remove a few others across phases.
                eligible = np.r_[np.arange(np_place), np.arange(np_place+np_object,n_cells)]
                present[rng.choice(eligible,max(1,n_cells//25),replace=False)] = False
                gids = np.flatnonzero(present); permuted = rng.permutation(gids)
                mapping[permuted,pi] = np.arange(len(permuted),dtype=float)
                fields['x'][pi] = xy[:,0]; fields['y'][pi] = xy[:,1]; fields['v'][pi] = v
                fields['t'][pi] = np.arange(n)/20; fields['S'][pi] = S[:,permuted]
            act = {k:fields[k] for k in ['t','S']}; act['mapping'] = mapping
            savemat(out/animal/(day+'.mat'), dict(pos={k:fields[k] for k in ['x','y','v']},
                     act=act, sess={'day_name':day}), do_compression=True)
            loaded = load_session(animal,day+'.mat',str(out))
            assert loaded['type'] == 'merged' and loaded['mapping'].shape == (n_cells,4)
            for pi, sub in enumerate(loaded['subsessions']):
                assert sub['S'].shape == fields['S'][pi].shape and len(sub['x']) == int(duration*20)
                ids = np.flatnonzero(np.isfinite(mapping[:,pi]))
                cols = mapping_columns(loaded['mapping'], ids, pi, sub['S'].shape[1])
                assert np.array_equal(sub['S'][:,cols], fields['S'][pi][:,mapping[ids,pi].astype(int)])
            print('Generada y verificada: %s/%s' % (animal,day),flush=True)
    pd.DataFrame(truth).to_csv(out/'ground_truth.csv',index=False)


def main():
    p = argparse.ArgumentParser(description='Datos sintéticos CA1')
    p.add_argument('--out-dir', required=True); p.add_argument('--seed',type=int,default=0)
    p.add_argument('--duration-s',type=float,default=240); p.add_argument('--n-cells',type=int,default=80)
    p.add_argument('--n-animals',type=int,default=2)
    a = p.parse_args(); generate(a.out_dir,a.seed,a.duration_s,a.n_cells,a.n_animals)


if __name__ == '__main__': main()
