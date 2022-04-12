# Using FDMT to Find a Radio Transient

## Background

In radio astronomy, radio transients (like fast radio bursts and pulsars) are a new and exciting field of study. Instead of imaging a supposed unchanging sky, radio astronomers are now looking for short timescale (from a few milliseconds to several seconds) signatures. For some of these phenomenon, like FRBs, the high-energy process that generates these signatures is not yet known. Detecting these signatures and combining with data from other telescopes and scientific efforts (like LIGO) may unlock brand new understandings to the nature of our universe.

For most of these radio transients, the source star is usually galaxies away. As the source generated a broadband RF pulse, the propagation of that pulse through ionized interstellar media disperses that pulse in time. This is commonly reffered to as the "sad trombone" effect as instead of recieving all the frequencies all at once, you recieve a descending sweep.

The shape of this sweep, or how quickly it descends is characterized by the dispersion measure (DM) of the pulse. Low DM pulses can have very little sweep, approaching a normal wideband pulse while high DM pulses can take seconds to sweep across all the bands. As the amount of dispersion is directly related to how much interstellar media it travelled through, we can conclude that the high DM pulses come from very far away.

It's easy to look at pulses once the dispersion has been removed ("dedispersed"), but it's difficult to find dispersed pulses in dynamic spectra. That's where libraries such as this one come in to play. By looking at many possible DMs over the chunk of dynamic spectra, we can find the one that maximizes the frequency-integrated power. That is to say, when there exists a dispersed pulse in the spectra, applying the correct DM transfromation to dedisperse and integrating across frequency will result in a power maximum.

This is what we want to find.

## Getting the Data

We'll be playing with data from Chris Bochenek's [STARE2](https://arxiv.org/abs/2001.05077) project. Specifically, from his fantastic Nature paper on [FRB 200428](https://www.nature.com/articles/s41586-020-2872-x).

First, we'll grab the published filterbank data from the Caltech library [here](https://data.caltech.edu/records/1647)
```julia
fb_file = download("https://data.caltech.edu/tindfiles/serve/8dcfd119-e2bc-458f-b56e-cb6ba10fd11f/")
```

Now we can use [SIGPROC.jl](https://github.com/kiranshila/SIGPROC.jl) to read out the data
```julia
using SIGPROC
fb = Filterbank(fb_file)
```

Might as well take a look at the raw dynamic spectra. Do you see a pulse in here? I certainly don't.
```julia
using Plots
plot(fb)
```

