"""Blocked cross-validated Poisson position decoder (Zhang et al., 1998).

Histograms smoothed in numerator and occupancy are the binned equivalent of
lab continuous Gaussian maps. Thresholds and cell selection use training only.
"""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import gaussian_filter
from .common import activity, shift_choices


def bin_phase(phase, width):
    t = phase['t']
    labels = np.floor((t - t[0] + 1e-9) / width).astype(int)
    n = int(np.floor((t[-1] + phase['dt'][-1] - t[0] + 1e-8) / width))
    good = labels < n
    aggregate = sparse.csr_matrix((np.ones(good.sum()), (labels[good], np.flatnonzero(good))),
                                   shape=(n, len(t)))
    nf = np.asarray(aggregate.sum(axis=1)).ravel()
    if np.any(nf == 0):
        raise ValueError('bins temporales vacíos; revisar time-bin/timestamps')
    xy = np.column_stack([aggregate @ np.nan_to_num(phase[k]) / nf for k in ['x', 'y']])
    speed = aggregate @ np.nan_to_num(phase['speed']) / nf
    valid = np.asarray(aggregate @ (~phase['valid']).astype(int)).ravel() == 0
    durations = np.asarray(aggregate @ phase['dt']).ravel()
    valid &= np.abs(durations - width) < max(.01, .1 * width)
    return dict(xy=xy, speed=speed, valid=valid, aggregate=aggregate, labels=labels,
                width=width, times=t[0] + (np.arange(n)+.5)*width, S=phase['S'])


def prepare_folds(binned, opts, eligible):
    n = len(binned['xy'])
    edges = np.linspace(-25, 25, int(np.ceil(50/opts.decoder_bin))+1)
    spacing = edges[1]-edges[0]
    centres = (edges[:-1]+edges[1:])/2
    xx, yy = np.meshgrid(centres, centres, indexing='ij')
    arena = xx*xx+yy*yy <= 625
    pixels = np.column_stack((xx[arena], yy[arena]))
    xy = binned['xy']
    ix = np.clip(np.searchsorted(edges, xy[:, 0], side='right')-1, 0, len(centres)-1)
    iy = np.clip(np.searchsorted(edges, xy[:, 1], side='right')-1, 0, len(centres)-1)
    folds = []
    for block in np.array_split(np.arange(n), opts.k_folds):
        if not len(block):
            continue
        test = block[eligible[block]]
        times = binned['times']
        train = eligible & ((times < times[block[0]]-opts.fold_gap_s) |
                            (times > times[block[-1]]+opts.fold_gap_s))
        if len(test) == 0:
            continue
        if train.sum() < 5:
            raise ValueError('insuficientes bins de entrenamiento')
        tr = np.flatnonzero(train)
        hist = sparse.csr_matrix((np.ones(len(tr)), (ix[tr]*len(centres)+iy[tr], tr)),
                                 shape=(len(centres)**2, n))
        occ_raw = np.asarray(hist @ np.full(n, opts.time_bin)).reshape(xx.shape)
        occ = gaussian_filter(occ_raw, opts.sigma/spacing, mode='constant')
        prior = np.ones(len(pixels)) if opts.prior == 'uniform' else occ[arena].copy()
        prior = np.maximum(prior, 1e-12); prior /= prior.sum()
        chance = np.ones(len(pixels)) if opts.prior == 'uniform' else occ_raw[arena].copy()
        chance /= chance.sum()
        folds.append(dict(train=train, test=test, hist=hist, occ=occ, prior=prior,
                          chance=chance, arena=arena, shape=xx.shape, sigma=opts.sigma/spacing))
    if not folds:
        raise ValueError('sin folds evaluables')
    return folds, pixels


def fold_counts(binned, train, opts, shift):
    # Shift the complete bin timeline, including pauses; do not concatenate running bins.
    n = len(train)
    source_bins = np.roll(train, -shift)
    frame_mask = (binned['labels'] < n)
    frame_mask &= source_bins[np.minimum(binned['labels'], n-1)]
    S = binned['S']
    values = activity(S, opts.activity, opts.threshold_sd, source=S[frame_mask])
    return np.roll(np.asarray(binned['aggregate'] @ values), shift, axis=0)


def cross_validate(binned, folds, pixels, opts, shift=0, rng=None):
    rows = []
    for fold_number, fold in enumerate(folds):
        counts = fold_counts(binned, fold['train'], opts, shift)
        keep = counts[fold['train']].sum(axis=0) >= opts.decode_min_events
        if not keep.any():
            raise ValueError('sin células activas en entrenamiento')
        counts = counts[:, keep]
        maps = np.asarray(fold['hist'] @ counts).T.reshape((-1,) + fold['shape'])
        smooth = gaussian_filter(maps, (0, fold['sigma'], fold['sigma']), mode='constant')
        rates = np.maximum(smooth[:, fold['arena']] /
                           np.maximum(fold['occ'][fold['arena']], 1e-12), 1e-6)
        test = fold['test']
        logp = counts[test] @ np.log(rates) - opts.time_bin*rates.sum(axis=0)[None, :]
        logp += np.log(fold['prior'])[None, :]
        prediction = pixels[logp.argmax(axis=1)]
        if shift != 0:
            rows.append(np.linalg.norm(prediction-binned['xy'][test], axis=1))
            continue
        posterior = np.exp(logp-logp.max(axis=1, keepdims=True))
        posterior /= posterior.sum(axis=1, keepdims=True)
        mean_prediction = posterior @ pixels
        chance = pixels[rng.choice(len(pixels), len(test), p=fold['chance'])]
        for j, idx in enumerate(test):
            rows.append(dict(bin_index=int(idx), time_s=binned['times'][idx], fold=fold_number,
                        true_x=binned['xy'][idx, 0], true_y=binned['xy'][idx, 1],
                        decoded_x=prediction[j, 0], decoded_y=prediction[j, 1],
                        speed_cm_s=binned['speed'][idx], n_cells_fold=int(keep.sum()),
                        error_cm=float(np.linalg.norm(prediction[j]-binned['xy'][idx])),
                        posterior_mean_error_cm=float(np.linalg.norm(mean_prediction[j]-binned['xy'][idx])),
                        chance_error_cm=float(np.linalg.norm(chance[j]-binned['xy'][idx]))))
    return np.concatenate(rows) if shift != 0 else pd.DataFrame(rows)


def decode_phase(binned, opts, rng, eligible=None):
    if eligible is None:
        eligible = binned['valid'] & (binned['speed'] >= opts.decode_min_speed)
    folds, pixels = prepare_folds(binned, opts, eligible)
    prediction = cross_validate(binned, folds, pixels, opts, rng=rng)
    shifts = shift_choices(len(binned['xy']), opts.time_bin, rng, opts.decode_shuffles)
    null = np.array([np.median(cross_validate(binned, folds, pixels, opts, shift=int(s))) for s in shifts])
    obs = float(prediction.error_cm.median())
    metrics = dict(n_cells=binned['S'].shape[1], n_cells_min_fold=int(prediction.n_cells_fold.min()),
                   n_bins=len(prediction), running_time_s=len(prediction)*opts.time_bin,
                   time_bin_s=opts.time_bin, min_speed=opts.decode_min_speed,
                   activity_mode=opts.activity, median_error_cm=obs,
                   mean_error_cm=prediction.error_cm.mean(),
                   posterior_mean_error_cm=prediction.posterior_mean_error_cm.mean(),
                   shuffle_mean_of_medians_cm=null.mean() if len(null) else np.nan,
                   shuffle_median_of_medians_cm=np.median(null) if len(null) else np.nan,
                   shuffle_p5_cm=np.percentile(null, 5) if len(null) else np.nan,
                   shuffle_p=(1+np.sum(null <= obs))/(len(null)+1) if len(null) else np.nan,
                   shuffle_z=(obs-null.mean())/null.std(ddof=1) if len(null)>1 and null.std()>0 else np.nan,
                   chance_median_error_cm=prediction.chance_error_cm.median())
    return metrics, prediction, null
