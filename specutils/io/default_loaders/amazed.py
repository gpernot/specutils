"""
Loader for Amazed spectrum files
"""
import os
import re

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.units import Unit, def_unit
from astropy.nddata import StdDevUncertainty

import numpy as np

from specutils.io.registers import data_loader, custom_writer
from specutils import Spectrum1D

__all__ = ['spec_identify', 'spec_loader']

_spec_pattern = re.compile(r'(?P<prefix>.*)_F\.fits')


def spec_identify(origin, *args, **kwargs):
    """
    Check whether given filename is FITS. This is used for Astropy I/O
    Registry.
    """
    return (isinstance(args[0], str) and
            _spec_pattern.match(args[0]) is not None)


@data_loader(label="Amazed spectrum", identifier=spec_identify, extensions=['fits'])
def spec_loader(file_name, **kwargs):
    """
    Loader for Amazed spectrum files.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]
    m = _spec_pattern.match(file_name)

    hdulist = fits.open(file_name, **kwargs)
    header = hdulist[0].header
    meta = {'header': header}

    # spectrum is in HDU 1
    data = hdulist[1].data['Flux']
    unit = Unit('1e-17 erg / (Angstrom cm2 s)')

    try:
        # Load companion error FITS file, if any
        error_file_name = '{}_ErrF.fits'.format(m['prefix'])
        errorhdu = fits.open(error_file_name, **kwargs)
    except:
        uncertainty = None
    else:
        error = errorhdu[1].data['Err']
        uncertainty = StdDevUncertainty(error)
        errorhdu.close()

    wave = hdulist[1].data['Wave']
    wave_unit = Unit('Angstrom')

    #mask = hdulist[1].data['and_mask'] != 0
    hdulist.close()

    return Spectrum1D(flux=data * unit,
                      spectral_axis=wave * wave_unit,
                      uncertainty=uncertainty,
                      meta=meta)
