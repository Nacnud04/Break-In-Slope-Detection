# Ice Plains

Simply put, ice plains are regions of lightly grounded ice between the ice sheet and ice shelf and are recognizable via their roughly constant downward sloping profile.

In this instance, the contact between the ice sheet and ice plain is known as the coupling line, and there is no official point for break in slope. 

For this use case, the same method is used to predict the location of the coupling line as the method to predict the break in slope. The two contacts can be defined by nearly the same parameters, so no modifications to the methodology need to be done. Additionally, it is assumed that the coupling line and the break in slope are contiguous.

Both the contact between the ice sheet and ice plain, as well as the ice plain and ice shelf display very similar patterns. This makes it extremely hard to differentiate between the two contacts, so ~25% of the time the plain-shelf contact is predicted as the coupling line.

## Detecting the ice plain from scratch

Accurately predicting the ice plain is a very significant challenge as the break in slope resembles the ice shelf, with a bunch of seemingly random spikes. Hence, to the algorithm it appears as an ice shelf with many break in slopes, so lots of possible so lots of possible Ib's get picked. This is used to the algorithms advantage for ice plain prediction.

To predict an ice plain, a span of data from the chosen break in slope, to the nearest other (non chosen) break in slope are selected. Then it checks to see if the average of the flowslope is just below zero, and the lower bound of the 2nd standard deviation from the mean is above a certain value. If these parameters are met, then the selection is expanded to the next closest (non chosen) break in slope. When the selection gets to its maximum possible size where the conditions are still true, then the algorithm is complete.

This methology is highly flawed. While it does correctly detect ice plains, it also detects a bunch of other data which is completely incorrect. This needs to be developed much more.

## Detecting the ice plain using pre-existing data

A quite accurate method of detecting the ice plain, and much more computationally efficent involves using the moa 2009 grounding line dataset. This dataset is based on more seaward grounding line features than the break in slope (specifically H- Hydrostatic equillibrium, and F- Point of tidal flexure).

As the methodology presented here predicts the most landward points (break in slope & coupling line), it is assumed that if there is significant distance between the coupling line and the moa 2009 dataset, and the moa dataset is seaward of the coupling line, then there is an ice plain.

This method works quite well for the area heavily studied here (Bungenstockrucken), however due to the age of the other dataset and lack of testing in other regions, this technique - while accurate in this case - is unreliable.
