mazer-jackson
=============

### Summary

Code base for V4 cell analysis written by Brandon Jackson for Jamie Mazer's
 lab

### Examples

Examples are contained in the `Examples/` directory. Each example consists
of a matlab file containing demo code and png files that show the output of
the sample code.

### Code Overview

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

**Classes**:

- AnglePlayFitter (@todo)
- GratRevFitter

### Loaders

Loaders are used to load stimulus images. All loaders have consistent
public APIs to make it easy for Models or other classes to interact with
multiple stimulus classes. All loaders must be passed a p2m file in the 
constructor so they can generate the correct stimuli.

For example, here is the code to load random stimuli for both an AnglePlay
stimulus and a GratRev stimulus.

    angleplay_pf = pffind('romeo0295*curvplay');
    APL = AnglePlayLoader(angleplay_pf);
    angleplay_stimulus = APL.randomStimulus();
    
    gratrev_pf = pffind('romeo0295*gratrev*001');
    GRL = GratRevLoader(gratrev_pf);
    gratrev_stimulus = GRL.randomStimulus();

These functions both use the public function `randomStimulus()`.

The stimuli generated by the loaders are downsampled versions of the ones
presented in the experiments. The amount of downsampling is controlled via
the `SuperLoader.DOWNSAMPLE_SCALE` constant. The stimulus image appears
centered in the frame. The stimuli are double matrices with range [0,1].

**Shared Public Methods**: 

- `randomStimulus()`
- `plotRandomStimuli()`

**Classes**:

- AnglePlayLoader
- GratRevLoader
- GratRevBarLoader
- GridCurvLoader
- SuperLoader
- WaveRev2Loader (@todo)

### Models

Models are used to generate responses to stimuli. All loaders have consistent
public APIs to make it easy for Writers or other classes to interact with
multiple models.

**Shared Public Methods**:

- `stimulate()`

**Classes**:

- ComplexCellModel
- GaborModel
- SimpleCellModel
- SuperModel

### Writers

Writers take p2m files and re-run experiments, presenting the stimuli
(loaded via a Loader class) to a Model, and then recording the responses
to a new synthetically-generated p2m file.

In addition to using a model to generate a spike train, it is also possible
to generate a poisson spike train using the `homogeneousPoisson()` and 
`inhomogeneousPoisson()` methods.

**Shared Public Methods**:

- `writeFromModel()`
- `discardRandomSpikes()`
- `homogeneousPoisson()`
- `inhomogeneousPoisson()`

**Classes**:

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

Finally there is a class that handles creating Gabor filters:

- GaborUtil

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

To-Do List
----------

- Convert batch scripts to load lists of files using the database, instead
of scanning directories filled with data files
- Clean up the code used to generate batch pstrf_delta results, and add to
 the repo
- I need to ensure that the ouputs of all model classes are consistent, and
then document it, perhaps in a how-to-write-a-model guide.
- A new system for organizing analysis code is necessary, because at the
moment it is very difficult to figure out whether code is in `RasterUtil`,
`SuperUtil`, or `PFUtil`. For example, perhaps `SuperUtil.autocorrelogram()`
should be moved to a class dedicated to correlograms? Should they be base
classes that the Core classes extend? For example right now I have defined 
a custom explainableVariance() method in AnglePlay but it should 
intelligently override the RasterUtil version somehow.
- `AnglePlayUtil` and `GratRevUtil` should be made static methods of their
respective Core classes
- `SuperUtil.bootstrappedAutocorrelogram()` should be merged into 
`SuperUtil.autocorrelogram()`, and perhaps moved to a new class
- Poisson analysis code should be moved out of `SuperUtil` and into its own
 class