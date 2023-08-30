# Break in slope detection and filtering algorithm

## 0. Removal of bad tracks & inadequate data
Initially for each track the necessary data is extracted, then the along track distance is computed. Next, the exact intersection location of the track with the grounding line is processed and along track distance is referenced relevant to the grounding line (0 = intersection with gline).

Next, distance to the grounding line is computed for each point, and a butterworth filter of order `5` and cutoff `0.032` is applied to the flowslope. This seems to be standard for applications such as this.

Tracks with an insufficent amount of points in the study area are removed.

Data is then interpolated over, and gaps in data (which were subsequently filled in by interpolation) of greater than `40m` are identified and removed again. Additionally, areas in the data where a value is repeated over and over for multiple measurements are removed, as this data is incorrect. This is done via computing standard deviation (less than `1e-10`, data is removed).

Tracks with more than `1/3` nan data are removed.

## 1. Initial location of slope break options
First, the derivative of the flowslope is calculated. This is referred to as the slope break.

The track is then split into segments which each contain 1 intersection with the grounding line model. This should allow for multiple intersections with the grounding line 

Possible slope break options are then chosen by looking for peaks in the slope break of a magnitude of `0.003` or greater. If less than 3 peaks are found, (common in areas with more gradual slopes) the magnitude is gradually reduced to `0.0001` until an adequate number of peaks are found.

Additionally, if multiple peaks are within `400` indicies of one another (about 3km), the peak with the highest magnitude is chosen. This ensures that usually only one point is selected within the grounding zone.

## 2. Capture parameters local to each peak for filtering
In order to filter out bad peaks from good peaks we capture parameters local to each peak to gauge its likelyhood of being the slope break.

First the percentage of Nan is computed for the points local to the selected peak. If this passes a threshold (currently set at `5%`) then the peak is omitted.

Additionally, points local to bad or "rough" data are removed. This is performed by removing points where more than `2.5%` of the local values fall outside of 2 standard deviations from the mean.

Then the following parameters are computed:
* $x$ Peak along track location
* $f$ Flowslope at peak
* $\sigma_0$ Standard deviation of flowslope local to peak (about `0.5km` radius)
* $\overline{f}_r$  $\overline{f}_l$ The mean of flowslope up and down track out to a certain radius (about `5km`)
* $\sigma_r$  $\sigma_l$ The standard deviation of flowslope both up and down track out to a certain radius (`40km`) while excluding points from an inner radius around the peak (`3km`)
* $qs$ The "quality score"

Quality score is a parameter to gauge the likelyhood of a peak being the actual break in slope. <br>
Mathmatically it is represented as:<br>
$$qs = \log{\frac{\sigma_r}{\sigma_l}} + \log{\frac{\overline{f}_r}{\overline{f}_l}} + 200\sigma_0 - 0.01x$$

<br>
**Deviation of heights:**
There is an option to compute the standard deviation of heights across all cycles on either side of the peak. Then the difference in the magnitude between the left and the right hand side will be added to the quality score, encouraging a large difference in deviation. This is done, as an ice shelf will have a larger deviation due to tidal influence than the sheet.

This system is currently flawed however, as it is extremely computationally expensive (need to open many other cycles worth of data), and the benefit is currently negligible.

## 3. Filter points by captured parameters
This is the stage where a point actually gets selected from the parameters which were just created.

Initially all peaks with a quality score of less than 2.5 are removed. Then filtering is applied to each point, which attempts to match certain properties. Specifically it ensures that one side of the peak has an average flowslope between `-0.001` and `0.002` (corresponding to the properties of an ice shelf) and the other side has an average flowslope less than `-0.002` (corresponding to a ice sheet).

The final peak is chosen as the point remaining with the higest quality score.
