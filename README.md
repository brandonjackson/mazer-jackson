mazer-jackson
=============

Code base for V4 cell analysis written by Brandon Jackson for Jamie Mazer's
 lab

There are two general kinds of classes in the code base: those that are 
objects designed to facilitate model building and testing and those that
 perform analysis. They are explained in turn below.

Model Building Classes
----------------------

### Philosophy

These classes contain several different kinds of objects with the 
intention of making it easier to perform analysis tasks like loading 
stimulus images, building models, fitting models, re-running an experiment 
with a model, and saving a model's spike train to a new p2m file.

### Fitters

Fitters are used to fit Models to a given p2m file.

- AnglePlayFitter (@todo)
- GratRevFitter

### Loaders

Loaders are used to load stimulus images. All loaders have consistent
public APIs to make it easy for Models or other classes to interact with
multiple stimulus classes.

- AnglePlayLoader
- GratRevLoader
- GratRevBarLoader
- GridCurvLoader
- SuperLoader (@todo)
- WaveRev2Loader (@todo)

### Models

Models are used to generate responses to stimuli. All loaders have consistent
public APIs to make it easy for Writers or other classes to interact with
multiple models.

- ComplexCellModel
- GaborModel
- SimpleCellModel
- SuperModel

### Writers

Writers take p2m files and re-run experiments, presenting the stimuli
(loaded via a Loader class) to a Model, and then recording the responses
to a new synthetically-generated p2m file.

    - AnglePlayWriter
    - GratRevWriter
    - SuperWriter

Analysis Classes
----------------

### Cores

These classes contain the bulk of the analysis code for each experiment 
type.

- AnglePlay
- GratRev
- GridCurv
- WaveRev2 (@todo)

### Utils

These classes perform analysis tasks that are shared (in theory) across 
multiple tasks.

    - PFUtil
    - RasterUtil
    - SuperUtil

There are two exceptions, which contain util functions specific to a
stimulus class. These might be merged into their respective Core classes
as static methods in the near future:

- AnglePlayUtil
- GratRevUtil

### Misc

- RatesDB: performs analysis of many cells and stores results in a searchable database format
- STRF: Spatio-Temporal Receptive Field analysis class

Scripts
-------

- Dependencies
    - sdf
- Interfaces
    - pangleplay
    - pgridcurv2
    - preel
    - pstrf
    - pstrf_delta
- Utils
    - pffind
    - gabor_filter
    - gabor_filter_translate

To-Do List
----------

- Loader classes should output consistent image sizes (e.g. 400 x 400) with 
consistent range (e.g. [0,1] or [-0.5,0.5])
- `gabor_filter.m` and gabor_filter_translate.m` should be moved to a class
definition
- `SuperUtil.bootstrappedAutocorrelogram()` should be merged into 
`SuperUtil.autocorrelogram()`, and perhaps moved to a new class