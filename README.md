# Antonis_Analysis
Code to analyse Antonis' data


Code to take data from 8 tetrodes (32 electrodes) implanted into brain of rats, and see how the principal components change with learning and memory over 6 sessions of a task. Observe place fields. Observe if activity changes in a fragile X rat model.

The information from the tetrodes can be used to triangulate location and provide limits of spatial resolution.

The principal component analysis is performed in KlustaKwik, providing cluster-identified spiking data and metadata.

The activity can be defined in terms of time and space to allow place field analysis.

Desire to combine all the sessions - 6 sessions, each ~10 minutes - and smooth the aligned data. Observe correlation between spatial bins in terms of spiking information.

Later, the desire is to compare with the LFP, to see if there's a difference in modulation between the WT and FMRP rats.

It is necessary to pre-define the experimentally relevant periods, to allow processing of reduced data set, that can be processed by a normal computer.

Central question - how well do different brains update their ensembles in the face of new information? How flexible are the neural circuits? Is the flexibility affected by Fragile X.

There is a desire to produce a power spectrum from the LFP data, and observe how this changes with time, and with genotype. To do this, we would need to exclude time periods when the animal is immobile.

There is a question as to whether learning and flexibility in neuronal ensembles is influenced by particular phases of activity, for example, are cells phase locked to gamma waves of neuronal activity.

The current experimental observation is that KO animals don't detect novelty, and this may be related to waves of activity in the hippocampus. WT spiking is phase locked during novelty, allowing flexibility of existing ensembles for new learning and memory formation.
