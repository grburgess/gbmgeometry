---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.8.0
  kernelspec:
    display_name: Python3
    language: python
    name: Python3
---

<!-- #region -->
# GBM Geometry Demo


**gbmeometry** is a module with routines for handling GBM geometry. It performs a few tasks:
* creates an astropy coordinate frame for Fermi GBM given a quarternion and spacecraft position
* allows for coordinate transforms from Fermi frame to an astropy frame (J2000, etc.)
* plots the GBM NaI detectors at a given time for a given FOV
* determines if an astropy SkyCoord location is within a NaI's FOV
* creates interpolations over GBM quarternions and SC coordinates

<!-- #endregion -->

```python
%matplotlib notebook
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u
from gbmgeometry import *
from gbmgeometry.utils.package_utils import get_path_of_data_file


```

## Interpolating the spacecraft position
First let's create an interpolating object for a given TRIGDAT file (POSHIST files are also readable)


```python
interp = PositionInterpolator.from_trigdat(get_path_of_data_file("glg_trigdat_all_bn080916009_v02.fit"))
```

The quaternion and sc_pos functions can take as inputs the time since trigger

```python
#
print ("Quaternions")
print (interp.quaternion(0))
print (interp.quaternion(10))
print
print ("SC XYZ")
print (interp.sc_pos(0))
print (interp.sc_pos(10))




```

## Single GBM detector properties

One can look at a single detector which knows about it's orientation in the Fermi SC coordinates as well as where it is currently pointing on its optical axis in J2000. In fact, since the GBM  frame is part of the astropy coordinate family, it can be transformed into any system!







```python
interp.quaternion(1)
```

```python
na = NaIA(interp.quaternion(1))
print (na.center)
print (na.center_icrs) #J2000
print (na.center.galactic) # Galactic
 
print ("Changing in time")
na.set_quaternion(interp.quaternion(100))

print (na.center)
print (na.center_icrs) #J2000
print (na.center.galactic) # Galactic
 
```

#### We can also go back into the GBMFrame 

```python
center_j2000 = na.center_icrs
center_j2000
```

```python
center_j2000.transform_to(GBMFrame(**interp.quaternion_dict(100.)))
```

### Earth Centered Coordinates

The sc_pos are Earth centered coordinates (in km for trigdat and m for poshist) and can also be passed. It is a good idea to specify the units!




```python
na = NaIA(interp.quaternion(0),interp.sc_pos(0)*u.km)
na.get_center()
```

## Working with the GBM class

Ideally, we want to know about many detectors. The GBM class performs operations on all detectors for ease of use. It also has plotting capabilities.

```python
myGBM = GBM(interp.quaternion(0),sc_pos=interp.sc_pos(0)*u.km)
```

We can either plot the detectors with a field of view:

```python
myGBM.plot_detector_pointings(fov=10);
```

or just where the optical axis is pointing

```python
myGBM.plot_detector_pointings(c='y', alpha=1,s=10);
```

```python
myGBM.plot_detector_pointings(fov=10, facecolor='r', projection = "astro globe", center = SkyCoord(30, -30, unit='deg', frame="icrs"), show_earth=False);
```

```python

myGBM.get_centers()


```

```python
[x.icrs for x in myGBM.get_centers()]
```

## Source/Detector Separation
We can even look at the separation angles for the detectors and a source.

```python
grb = SkyCoord(ra=130.,dec=-45 ,frame='icrs', unit='deg')

seps = myGBM.get_separation(grb)

seps.sort("Separation")

seps
```

# Fermi plotting

```python
from gbmgeometry.spacecraft.fermi import *


f = Fermi(quaternion=[0,0,0,0], sc_pos=np.array([1,1,1]))
f.plot_fermi(color_dets_different=True, plot_det_label=False);
```

```python

```
