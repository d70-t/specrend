# -*- coding: utf-8 -*-
"""
Color rendering of spectra
==========================

:author: Tobias Kölling <tobias.koelling@physik.uni-muenchen.de>

This module can convert spectra to RGB. It is inspired by
John Walkers ``specrend.c`` and the respective python version
of Andrew Hutchins, but is completely rewritten.

The wavelenth resampling code is from Florian Ewald <florian.ewald@campus.lmu.de>.

Usage::

    >>> spectralimage = ... #some spectral image
    >>> spectralbands = ... #the corresponding wavelengths
    >>> converter = SpecRGBConverter(spectralbands)
    >>> rgbimage = converter.spectrum_to_rgb(spectralimage)

The original specrend.py comment::

                Colour Rendering of Spectra

                    by John Walker
                http://www.fourmilab.ch/

        Last updated: March 9, 2003

    Converted to Python by Andrew Hutchins, sometime in early
    2011.

        This program is in the public domain.
        The modifications are also public domain. (AH)

    For complete information about the techniques employed in
    this program, see the World-Wide Web document:

            http://www.fourmilab.ch/documents/specrend/

    The xyz_to_rgb() function, which was wrong in the original
    version of this program, was corrected by:

        Andrew J. S. Hamilton 21 May 1999
        Andrew.Hamilton@Colorado.EDU
        http://casa.colorado.edu/~ajsh/

    who also added the gamma correction facilities and
    modified constrain_rgb() to work by desaturating the
    colour by adding white.

    A program which uses these functions to plot CIE
    "tongue" diagrams called "ppmcie" is included in
    the Netpbm graphics toolkit:
        http://netpbm.sourceforge.net/
    (The program was called cietoppm in earlier
    versions of Netpbm.)

"""
import numpy as np

from .cietables import CIE1931


class ColorSystem(object):
    """
    Baseclass for a general color system.
    """
    name = 'UNDEFINED'
    gamma = None

    def add_gamma(self, rgb):
        """
        Adds gamma correction to a given rgb image.

        :note: The image must be given in float, normalized to 0...1
        """
        if callable(self.gamma):
            return self.gamma(rgb)  # pylint: disable=not-callable
        else:
            return rgb**(1.0 / self.gamma)

    def __repr__(self):
        return "<{} ColorSystem>".format(self.name)


class SpecrendColorSystem(ColorSystem):
    """
    ColorSystem corresponding to the original specrend implementation
    of John Walker.
    """
    def __init__(self, csdef):
        self.xyz2rgb = self.cs2mat(csdef)
        self.gamma = csdef["gamma"]
        self.name = csdef["name"]

    @classmethod
    def cs2mat(cls, cs):
        """
        Given an additive tricolour system CS, defined by the CIE x
        and y chromaticities of its three primaries (z is derived
        trivially as 1-(x+y)), and a desired chromaticity (XC, YC,
        ZC) in CIE space, determine the contribution of each
        primary in a linear combination which sums to the desired
        chromaticity.  If the  requested chromaticity falls outside
        the Maxwell  triangle (colour gamut) formed by the three
        primaries, one of the r, g, or b weights will be negative.

        Caller can use constrain_rgb() to desaturate an
        outside-gamut colour to the closest representation within
        the available gamut and/or norm_rgb to normalise the RGB
        components so the largest nonzero component has value 1.
        """
        # pylint: disable=invalid-name,too-many-locals
        xr = cs["xRed"]
        yr = cs["yRed"]
        zr = 1 - (xr + yr)
        xg = cs["xGreen"]
        yg = cs["yGreen"]
        zg = 1 - (xg + yg)
        xb = cs["xBlue"]
        yb = cs["yBlue"]
        zb = 1 - (xb + yb)
        xw = cs["xWhite"]
        yw = cs["yWhite"]
        zw = 1 - (xw + yw)

        rx = (yg * zb) - (yb * zg)
        ry = (xb * zg) - (xg * zb)
        rz = (xg * yb) - (xb * yg)
        gx = (yb * zr) - (yr * zb)
        gy = (xr * zb) - (xb * zr)
        gz = (xb * yr) - (xr * yb)
        bx = (yr * zg) - (yg * zr)
        by = (xg * zr) - (xr * zg)
        bz = (xr * yg) - (xg * yr)

        rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw
        gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw
        bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw

        rx /= rw
        ry /= rw
        rz /= rw
        gx /= gw
        gy /= gw
        gz /= gw
        bx /= bw
        by /= bw
        bz /= bw

        return np.array(((rx, ry, rz), (gx, gy, gz), (bx, by, bz)))


class MatrixColorSystem(ColorSystem):
    """
    ColorSystem defined by a linear transformation w.r.t. XYZ space.
    """
    def __init__(self, name, xyz2rgb, gamma):
        self.name = name
        self.xyz2rgb = xyz2rgb
        self.gamma = gamma


def gamma_rec709(rgb):
    """
    Compute REC709 Gamma (according to the original specrend.py)
    """
    cc = 0.018
    return np.where(rgb < cc,
                    rgb * (((1.099 * (cc**0.45)) - 0.099) / cc),
                    (1.099 * (rgb**0.45)) - 0.099)


def gamma_srgb(rgb):
    """
    Compute sRGB Gamma according to:
    http://www.w3.org/Graphics/Color/sRGB

    :param rgb: rgb-array WITHOUT gamma, values 0.0...1.0 (linear data)

    :return: rgb-array WITH gamma, values 0.0...1.0 (nonlinear data)
    """
    return np.where(rgb <= 0.00304, 12.92*rgb, 1.055*(rgb**(1./2.4))-0.055)


def gamma_srgb_reverse(rgb):
    """
    Compute sRGB Gamma according to:
    http://www.w3.org/Graphics/Color/sRGB

    :param rgb: rgb-array WITH gamma, values 0.0...1.0 (nonlinear data)

    :return: rgb-array WITHOUT gamma, values 0.0...1.0 (linear data)
    """
    return np.where(rgb <= 0.03928, rgb/12.92, ((rgb+0.055)/1.055)**2.4)


# sRGB table is also from: http://www.w3.org/Graphics/Color/sRGB
CS_SRGB = MatrixColorSystem("sRGB",
                            np.array(((3.2410, -1.5374, -0.4986),
                                      (-0.9692, 1.8760, 0.0416),
                                      (0.0556, -0.2040, 1.0570))),
                            gamma_srgb)

# the following color systems stem from the original specrend.py
CS_NTSC = SpecrendColorSystem({"name": "NTSC",
                               "xRed": 0.67, "yRed": 0.33,
                               "xGreen": 0.21, "yGreen": 0.71,
                               "xBlue": 0.14, "yBlue": 0.08,
                               "xWhite": 0.3101, "yWhite": 0.3163,
                               "gamma": gamma_rec709})

CS_EBU = SpecrendColorSystem({"name": "SUBU (PAL/SECAM)",
                              "xRed": 0.64, "yRed": 0.33,
                              "xGreen": 0.29, "yGreen": 0.60,
                              "xBlue": 0.15, "yBlue": 0.06,
                              "xWhite": 0.3127, "yWhite": 0.3291,
                              "gamma": gamma_rec709})

CS_SMPTE = SpecrendColorSystem({"name": "SMPTE",
                                "xRed": 0.63, "yRed": 0.34,
                                "xGreen": 0.31, "yGreen": 0.595,
                                "xBlue": 0.155, "yBlue": 0.07,
                                "xWhite": 0.3127, "yWhite": 0.3291,
                                "gamma": gamma_rec709})

CS_HDTV = SpecrendColorSystem({"name": "HDTV",
                               "xRed": 0.67, "yRed": 0.33,
                               "xGreen": 0.21, "yGreen": 0.71,
                               "xBlue": 0.15, "yBlue": 0.06,
                               "xWhite": 0.3127, "yWhite": 0.3291,
                               "gamma": gamma_rec709})

CS_CIE = SpecrendColorSystem({"name": "CIE",
                              "xRed": 0.7355, "yRed": 0.2645,
                              "xGreen": 0.2658, "yGreen": 0.7243,
                              "xBlue": 0.1669, "yBlue": 0.0085,
                              "xWhite": 0.3333333333, "yWhite": 0.3333333333,
                              "gamma": gamma_rec709})

CS_REC709 = SpecrendColorSystem({"name": "CIE REC709",
                                 "xRed": 0.64, "yRed": 0.33,
                                 "xGreen": 0.30, "yGreen": 0.60,
                                 "xBlue": 0.15, "yBlue": 0.06,
                                 "xWhite": 0.3127, "yWhite": 0.3291,
                                 "gamma": gamma_rec709})


class CieColorTable(object):
    """
    Represents a color sensitivity function for given wavelenths.
    :param wvlns: list of wavelengths to interpolate the sensitivity on.
    """
    wvlns, table = map(np.array, zip(*CIE1931))

    def __init__(self, wvlns):
        compare = (wvlns != self.wvlns)
        try:
            compare = any(compare)
        except TypeError:
            pass
        if compare:
            self.resample_wvlns(wvlns)
        else:
            self.wvlns_idx = np.arange(0, len(self.wvlns), 1)
        self.nwvlns = len(self.wvlns)

    def resample_wvlns(self, wvlns):
        """
        Resamples the color table to a set of given wavelenths.

        :author: Florian Ewald <florian.ewald@campus.lmu.de>

        :param wvlns: Wavelengths to resample to.
        """

        self.wvlns_idx, = np.where((self.wvlns.min() < wvlns) & (wvlns < self.wvlns.max()))
        wvlns_common = wvlns[self.wvlns_idx]
        delta_wvlns_common = np.abs(wvlns_common[1:] - wvlns_common[:-1])
        delta_wvlns_common = np.hstack(([delta_wvlns_common[0], ], delta_wvlns_common))

        cie_colour_match = np.zeros((len(wvlns_common), 3), float)
        for i in range(3):
            cie_colour_match[:, i] = np.interp(wvlns_common, self.wvlns, self.table[:, i]) * delta_wvlns_common
            # ensure that the integrals of the colors stay equal
            cie_colour_match[:, i] /= cie_colour_match[:, i].sum()

        self.wvlns = wvlns_common
        self.table = cie_colour_match


def norm_array(arr):
    """
    Returns an array scaled to a maximum value of 1.
    """
    return arr/arr.max()


def array_max_dynamic_range(arr):
    """
    Returns an array scaled to a minimum value of 0 and a maximum value of 1.
    """
    finite_arr = arr[np.isfinite(arr)]
    low = np.nanmin(finite_arr)
    high = np.nanmax(finite_arr)
    return (arr - low)/(high - low)


def float2byte(arr):
    """
    Convert float [0 ... 1]-valued array to uint byte array.
    """
    return (arr*255.).astype('uint8')


def set_outliers_to_zero(arr):
    """
    Set all pixels which are bigger than 100 times the median to 0
    """
    arr[arr >= np.median(arr)*100.] = 0.


def constrain_rgb(rgb, inplace=False):
    """
    If the requested RGB shade contains a negative weight for
    one of the primaries, it lies outside the colour gamut
    accessible from the given triple of primaries.  Desaturate
    it by adding white, equal quantities of R, G, and B, enough
    to make RGB all positive.
    """
    white = np.nanmin(rgb)
    if white < 0:
        if inplace:
            rgb -= white
            return rgb
        else:
            return rgb - white
    else:
        return rgb


def postProcessRGB(rgb, colorsystem=CS_SRGB):
    """
    Post process an rgb image.

    This is normalize to full dynamic range, apply gamma and convert to 8bit uint.

    :note: This function will give errorneous results, if applied only on image slices.

    :note: This function operates in-place, so the input rgb will be overwritten!

    :param cs: Color system to use
    :param rgb: Preprocessed rgb image
    """
    constrain_rgb(rgb, inplace=True)
    set_outliers_to_zero(rgb)
    return float2byte(colorsystem.add_gamma(array_max_dynamic_range(rgb)))


class SpecRGBConverter(object):
    """
    Converter for spectra to RGB images.

    The default color system is sRGB because it is the most widely adopted
    color system on the web.

    :note: Check self.colortable.wvlns_idx for the taken spectral indices.
           The spectra will be contrained to these indices when the
           transformation is applied.

    :param wvlns: List of wavelengths of the spectra to convert.
    :param cs: RGB color system to convert to (default: SRGB)
    """
    def __init__(self, wvlns, colorsystem=CS_SRGB):
        self.colortable = CieColorTable(wvlns)
        self.colorsystem = colorsystem
        self.spec2rgb = np.tensordot(self.colortable.table,
                                     self.colorsystem.xyz2rgb,
                                     axes=((1,), (1,)))

    def crop_spectrum(self, spectrum):
        """
        Crop a given spectrum to the used wavelengths, assuming the spectrum
        is either given with the same wavelenths as in ``__init__`` or already
        cropped (-> this function is idempotent)

        :param spectrum: over-defined spectrum
        :return: spectrum only on needed wavelengths
        """
        nwvlns = spectrum.shape[-1]
        if nwvlns < self.colortable.nwvlns:
            raise ValueError(
                    'converter is configured for {} wavelenths, but I got only {}!'
                    .format(self.colortable.nwvlns, nwvlns))
        elif nwvlns > self.colortable.nwvlns:
            # assume shape of ``spectrum`` is according to ``wvlns`` in __init__ and cut.
            return spectrum[..., self.colortable.wvlns_idx]
        return spectrum

    def xyz_to_rgb(self, xyz):
        """
        convert xyz to RGB according to the defined color system
        :param xyz: Input image in xyz coordinates
        :return: Image in rgb coordinates (not normalized)
        """
        return np.tensordot(xyz, self.colorsystem.xyz2rgb, axes=((-1,), (1,)))

    def spectrum_to_xyz(self, spectrum):
        """
        Calculate the CIE X, Y, and Z coordinates corresponding to
        a light source with spectral distribution given by ``spectrum``.
        ``spectrum`` must be a ``numpy`` array and the last dimension
        must be the spectral dimension with equally spaced wavelengths
        equal to ``self.wvlns``.

        :param spectrum: Input spectrum
        :return: Image in xyz coordinates (not normalized)
        """

        # NOTE: Do NOT normalize XYZ, otherwise the gray information is lost!
        #       The image must be normalized globally after color conversion.
        return np.tensordot(self.crop_spectrum(spectrum), self.colortable.table, axes=((-1,), (0,)))

    def spectrum_to_rgb(self, source, postprocess=True):
        """
        Calculate the R, G, B coordinates corresponding to
        a light source with spectral distribution given by ``source``.
        The spectral axis must be the last axis, the spectrum
        must be defined in equidistant steps, either on the wavelenths
        defined during the ``__init__`` call, or already cropped to
        ``self.colortable.wvlns``.

        :note: Use ``postprocess`` only on full images, not on slices!

        :param source: a [...,nwvlns] - shaped input spectral distribution
        :param postprocess: If set to ``True``, a 'ready-to-store' rgb image
                            will be returned.
        :return: RGB image
        """
        rgb = np.tensordot(self.crop_spectrum(source), self.spec2rgb, axes=((-1,), (0,)))
        if postprocess:
            return postProcessRGB(rgb, self.colorsystem)
        else:
            return rgb
