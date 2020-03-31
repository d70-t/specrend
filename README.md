# Color rendering of spectra

This module can convert spectra to RGB. It is inspired by
John Walkers ``specrend.c`` and the respective python version
of Andrew Hutchins, but is completely rewritten.

Usage:

```python
spectralimage = ... #some spectral image
spectralbands = ... #the corresponding wavelengths
converter = SpecRGBConverter(spectralbands)
rgbimage = converter.spectrum_to_rgb(spectralimage)
```
