# A-New-Composite-Design-for-Chemometrics-and-Computer-Experiments
Related MATLAB codes for case studies and applications
1. Put the files into the same directory of MATLAB
2. Run models.m first to load all underlying models. Note that f4 is defined by f4.m.
3. Each case study of our manuscript is conducted by f#comparison.m file respectively. The procedure follows as:
    - Give initial designs as matrices copied from the existing literature and website.
    - Transform the design entries into practical experimental domain.
    - Derive experimental results of each initial design and each possible composite design and 1000 new random observations.
    - Derive Kriging models of each design and its MSE for predicting the random observations.
    - Codes with comment %quadratic can be ignored.
