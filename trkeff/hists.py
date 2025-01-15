from get_hists.x_h import *
from get_hists.y_h import *
from get_hists.xz_h import *
from get_hists.yz_h import *
from get_hists.chi2_h import *
from get_hists.chi2ndf_h import *
from get_hists.trkP_h import *
from get_hists.n_h import *
from get_hists.d0_h import *
from get_hists.e_h import *

xy_range = {
    'min': {'x': -70., 'y': 0. },
    'max': {'x':  10., 'y': 80.}
}
xy_eff_range = {
    'min': {'x': -42., 'y': 19.},
    'max': {'x': -10., 'y': 48.}
}

get_hists: dict = {
    'x':         get_h_x_eff,
    'y':         get_h_y_eff,
    'xz':        get_h_xz_eff,
    'yz':        get_h_yz_eff,
    'chi2':      get_h_chi2_eff,
    'chi2ndf':   get_h_chi2ndf_eff,
    'trkP':      get_h_trkP_eff,
    'n':         get_h_n_eff,
    'd0':        get_h_d0_eff,
    'e':         get_h_e_eff,

    'x.y':       get_h_x_y_eff,
    'x.xz':      get_h_x_xz_eff,
    'x.yz':      get_h_x_yz_eff,
    'x.chi2':    get_h_x_chi2_eff,
    'x.chi2ndf': get_h_x_chi2ndf_eff,
    'x.trkP':    get_h_x_trkP_eff,
    'x.n':       get_h_x_n_eff,

    'y.xz':      get_h_y_xz_eff,
    'y.yz':      get_h_y_yz_eff,
    'y.chi2':    get_h_y_chi2_eff,
    'y.chi2ndf': get_h_y_chi2ndf_eff,
    'y.trkP':    get_h_y_trkP_eff,
    'y.n':       get_h_y_n_eff,

    'xz.x':       get_h_xz_x_eff,
    'xz.y':       get_h_xz_y_eff,
    'xz.yz':      get_h_xz_yz_eff,
    'xz.chi2':    get_h_xz_chi2_eff,
    'xz.chi2ndf': get_h_xz_chi2ndf_eff,
    'xz.trkP':    get_h_xz_trkP_eff,
    'xz.n':       get_h_xz_n_eff,

    'yz.x':       get_h_yz_x_eff,
    'yz.y':       get_h_yz_y_eff,
    'yz.chi2':    get_h_yz_chi2_eff,
    'yz.chi2ndf': get_h_yz_chi2ndf_eff,
    'yz.trkP':    get_h_yz_trkP_eff,
    'yz.n':       get_h_yz_n_eff,

    'dxRef':      get_h_dxRef,
    'dyRef':      get_h_dyRef,

    'dxz':        get_h_dxz,
    'dyz':        get_h_dyz
}

prpts_data_all = (
    'x', 'y', 'xz', 'yz', 'chi2', 'chi2ndf', 'trkP', 'n', 'd0',
    'x.y', 'x.xz', 'x.yz', 'x.chi2', 'x.chi2ndf', 'x.trkP', 'x.n',
    'y.xz', 'y.yz', 'y.chi2', 'y.chi2ndf', 'y.trkP', 'y.n',
    'xz.x', 'xz.y', 'xz.yz', 'xz.chi2', 'xz.chi2ndf', 'xz.trkP', 'xz.n',
    'yz.x', 'yz.y', 'yz.chi2', 'yz.chi2ndf', 'yz.trkP', 'yz.n',
    'dxRef', 'dyRef', 'dxz', 'dyz'
)


def create_h(
    run_or_mcSet = 7080,
    track_types = (1, 11, 3, 13),
    prpts = prpts_data_all
):
    h = {}
    for tt in track_types:
        h[tt] = {}

        for prpt in prpts:
            h[tt][prpt] = get_hists[prpt](run_or_mcSet, tt)
    return h
