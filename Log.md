# to-do
* 

# done
* Post processing removes any detected point more than 10km from any other estimated points.
* Added variable break in slope threshold change.
  * Essentially the algorithm acts as normal, except it requires 3 or more peaks to be detected to choose the best one from. 
  * If this requirement is not met it decreases the threshold which counts as a peak from 0.003 to 0.002 to 0.001
* Attempted track segmentation based on grounding line intersection where the actual properties of each point are determined by the entire (not segmented) track.
  * This causes something weird to happen where little data is picked up.
* Added quality paramter. I tried making this parameter select the point entirely on it's own, but noisy and bad data is not filtered out. Therefore it sometimes happens to match and be of "good quality"
  * So I use this to supplement the previous procedure. Essentially I weakened the selection process and made it more forgiving, to allow for more points to pass through. To this I implemeted the quality parameter which makes the final selection on which point to use.
* I noticed a lot of bad data causes points to be way out in the ocean. This could just be natural surface roughness or something.
  * this causes a lot of points to be chosen as the possible break in slope however this is not the case.
  * To fix this I did the following:
    * Any points with more than 5% nan within the local area (defined by local radius, usually around 0.5-1km on either side) are omitted.
    * Any points with more than 5% of local values outside of 2 standard deviations from the mean are also removed.



## OLD STUFF
* Fix cheat on along track distance direction computation
  * Basically if the along track distance increases in one direction then we want to look for peaks with a certain sign
  * If the along track distance does the opposite we want to look for the opposite sign.
  * Right now I am purely figuring this out based on the polygon shape.
* ~~Do not perform Ib calculations based on a low angle between track and gline.~~
* ~~Reprocess all of the data so it is actually in order~~
* ~~Fix the ordering issue~~
  * ~~Thinking of reordering the data by the along track distance. But then why are the times so fucked?~~
* Look at flowslope on plot with deviations in color.
* ~~Add a parameter to notify if a track could be inaccurate due to a high angle between flow direction and track direction.~~
  * I don't believe this causes any issues whatsoever.
---
* ~~Add something which takes the average to the right, and keeps/removes a peak depending on it's value~~
  * ~~For instance, if the peak would normally be kept, but the average value of all points 2km to the right is greater than -0.002 then get rid of it.~~
* ~~Correct elevation average computation~~
* ~~Add 2nd derivative peak selection bounds.~~
  * ~~One should be put in place which is based around the offset from zero of the elevation of the data~~
    * ~~This can maybe be taken as the absolute value of the average being greater or less than a certain value?~~
  * One should be based off of the spikyness of the data
    * This does not seem to be needed yet, but this spikyness could definitely pose a threat.
* ~~Fix peak index to peak distance calculation
  * ~~The problem here is that it assumes that the track does not shift. So if the track does shift it assumes there are the same amount of points in the study area which there are not.~~
* ~~Remove tracks/cycles with less than 50% coverage. These are unreliable~~
* ~~Remove positive peaks from calculations~~
* ~~Add correction to do correlation allowing for nans~~
  * ~~To do this, I am filling in all nan locations with the linearlly interpolated, but not nan removed data.~~
  * ~~Then I substitute in the nan values~~
* ~~Add correction to allow stanard deviation calculations to work with tracks with missing points~~
* ~~Understand what is new in ATL06 v5 vs v6~~
* ~~Add correction to hide the following:
  * ~~/home/jovyan/GLine-Mig/ICESat2GroundingLineMigration/IceSatHDF5Unpacker.py:915: RuntimeWarning: invalid value encountered in double_scalars angle = math.atan(dy/dx)~~
* ~~When interpolate fill large gaps of missing data with nan. Aka do not interpolate over them~~
* ~~add correction for if grounding line is not found~~
* ~~Why on earth is data not being found for a certain cycle when I exported it previously??~~
  * ~~**THEORY:** The upgrade of the data from atl06 v5 to v6 could have caused this. Maybe not all the data has been processed and is on nasa earthdata yet for some reason?~~
  * ~~Previously I had exported data which returned many parameters within a polygon for a given date range. Eventually this stopped giving back any data. So for future ease and understandability I started exporting data by the cycle number. However the same exact period still doesnt return any data!!~~
  * ~~Additionally, during a cycle when there is data it returns nowhere near the amount of granules it should! Usually get around 65, but sometimes it only returns like 20!~~
  * ~~similar issues:~~
    * ~~does 2020-04-01 through 2020-06-30 really have no usable data? (cyele 11)~~
    * ~~does 2020-10-01 through 2020-12-21 have 0 usable granules in the study area? (cycle 13)~~
* ~~why is there multiple intersections when there definitely shouldn't be??~~
  * ~~This is because the gline is just a list of x and y points.~~
  * ~~Therefore if there is multiple lines then they are joined by straight lines and false intersections can be created~~
  * ~~To correct for this we need to create multiple linestrings.~~
* add ability to work with multiple intersections and only choose one
  * **What if we keep both, and plot the along track distance and use a vertical line to indicate a gline intersection**
  * Why is there so many multipoint results? Each time it nearly always returns a multipoint.
  * This has partially been done. However there needs to be a way to actually use the multiple intersections.
* ~~rewrite intersection script.~~
* ~~what to do about multiple grounding line intersections?~~
* ~~reexport data with datetime~~
* ~~numerically report how the cross correlation reduces the deviation~~
* ~~explain notebook~~
* ~~add cycle number to the output dataframe (parameter called `cycle_number` in group `orbit_info`)~~
* ~~exported data automatically gets split up by rgt and ground track~~
* ~~export more data over a longer time 2020~~
* ~~understang cloud = True~~ Damn I'm dumb
* ~~export data for 2022~~
* read to get some ideas on why the difference between the left and right track
* ~~capture deviation change
* ~~somehow get rgt 559 to output usable data???~~ **I think this one is impossible**
* ~~increase size of bounding polygon~~
* ~~re-export all data to get h_li~~
* apply standard deviation plot to rgt 559 (or nearby tracks)
  * Look at height along the plot, observe correlations between deviation height and flowslope

# Notable things done

* Plotted the assumed point F of tidal flexure against the offset of the elevation data from normal. (to represent tidal differences, somewhat similar to what is done in most papers)
  * My data aligns exactly with what is reported by [https://doi.org/10.5194/tc-2022-265](https://doi.org/10.5194/tc-2022-265)
  * The data even aligns with the tracks which dissapeared with the new version of ICESat2 ATL06 (v6)
* Noticed that the crappier tracks, (ones missing a lot of data), in sections with good data have much higher amplitude swings than the rest of the tracks.
  * Like the peaks are in the same place, just are of a much higher magnitude.
* Managed to get data from cycles 10, 11, 12, & 13. For some reason after the new update of the data version from v5 to v6 these dissapeared. 
  * Only 2 weeks later they are back however!
  * So i exported all of this 6/5/23
  * However they are not all back in full? Many versions of cycle 13 which were great are still MIA.
  * Also data is still much more crap in some places.
  * The data which is good is ever so slighly different than previous data. Maximum deviation appears to be around 0.0001 or less.
  * Considering the magnitude of the actual slope data this does not make any difference whatsoever.
* Added visualizations of individual point f's as well as peak recognition utilities
* Added removal of tracks which are lacking in data
  * So, when gaps of missing data are removed, a calcuation is performed to determine whether the track is worth keeping
  * If there is less than 66% of the entire track which remains then the entire track is removed.
* Added selective choosing of a peak on the 2nd derivative of heights
  * I think this could work well as a method to find point f?
  * By choosing where there is usually the lowest peak this most usually returns a decent value
  * Occassionally however things get tripped up. An example of this is 1062gt2r
  * Additionally some extra filtering was added which includes:
    * If there are multiple peaks within 100m of one another, the peak with the highest prominence is chosen
    * 95% of the data around a peak must not be NAN
    * Flowslope at the peak must be greater than (0.001?, 0.0015?, 0.002?)
  * Lastly this required the addition of a index to along track distance conversion.
* Fixed detection of multiple intersections
  * Since the grounding line was parsed into a list of x and y's if there was any gap in the grounding line where it is not in the study region a straight line was created straight through
  * This then creates multiple intersections
  * To fix this I added a calculation which splits the list of x and y into multiple line segments depending on if the gap between two datapoints is greater than a certain distance.
  * **This method will not work if there are truly multiple intersections with the grounding line**
* Corrected for certain rare edge cases in the flowslope calculations
  * This includes stuff such as when computing an angle of the track relative to polar stereographic north, rare angles such as pi / 2 result in a divide by zero error.
* Added correction to allow for computation of standard deviation for areas where only set amounts of data may be present.
  * This yields some interesting results, as the final deviation curve jumps between two or more different curves depending on the quantity of data it can compute the deviation for
* Added function which removes sections of unnaturally straight data.
  * Calculates standard deviation over segments.
* Threw higher level processing stuff (listed below) into a utility file
  * Includes: Nan removal, bad track removal, function which clumps grounding line intersection computations, cross correlations
* Added funct which removes large gaps in data (replaces with nan) which get filled in during a linear correlation.
  * Nan values make the correlation lag always the same no matter the data
    * Removal of nans does not change the correlation lag whatsoever (not even by 1m, it is the exact same). So it is quite possible that scipy does a linear interpolation on the data if there is a gap during the cross correlation process
  * This method did not result in the removal of tracks with a few points, but very few points in the study area. In this case the line in the study area is nearly perfectly straight. To correct for this I did the following:
    * Calculate the standard deviation of each track. For every single track (even ones with wildly different profiles) the standard deviation was around 0.0015 to 0.0027.
    * For the straight ish line track the deviation was on the magnitude of 10^-19.
    * This works as a great way to filter out these anomalies
* Noticed a strange pattern in rgt1062 gt3r & gt3l.
  * Cycle 8 & 17 are roughly the same, and the majority of the other tracks are also roughly the same. For gt3l cycles 8 & 17 take a longer time to experience grounding zone features along track. Where as for gt3r this takes a shorter along track distance to occur.
  * Maybe this indicates that the grounding line behaves quite differently just a few meters over, indicating something about subsurface features??
* Rewrote the part which works with the icepyx and nasa earthdata api.
  * This corrects for the dumb effect of exporting not working when cropped by cycle or date
  * Also improves program efficency somewhat.
  * **This still does not return the periods of data which appear to have just dissapeared??**
* Added correction for different kinds of intersections
  * These intersection types include: Point (1 intersection), Multipoint (Multiple intersections), Polyline (No intersections??)
* Added correction for if no intersection with the grounding line is found.
* Added extraction by cycle number
* Added QOL features to data extraction.
  * Allows for a better understanding of the quality of data being extracted
  * Allows for failures in data processing to occur, and the extraction continues
  * Notifies on RGT's with insufficient amounts of data.
* Exported data for 2020
* Analyzed plenty of other parameters such as time, xy position, etc. to get a better understanding of why one beam in the beam pair returns something different than the other beam.
  * *Status update:* Still have no fuckin clue man.
* Upgraded data inital processing program to only extract desired rgt's
* Implemented automatic data extraction into separate folders and files
* Cleaned up directory tree and renamed some notebooks
* Cleaned up rgt analysis notebook
  * Chucked a ton of code into a utility file
* Added documentation to the file with the most processing.
* Calculated the deviation over distance for both before and after the cross correlation is applied
  * Also simplified the computing of these by merging similar processes
* Computing of the along track distance for some reason removed the datetime parameter?
* Using along track distance, with 0 at the grounding line instead of distance from the grounding line
  * This is because the grounding line geometry near rgt 559, rgt 1062, and rgt 62 becomes parallel to the actual track
  * This causes discontinuities and major squishing in the data which causes all data processing to essentially fail.
* Discovered major slope variation in rgt 1062. And it cross correlates perfectly over a distance of 7km.
  * RGT 1062 has an extremely obvious grounding zone.
* Noted a major difference between the slope profile between a left and right beam. 
  * So I read a bunch and found nothing.
  * Eventually tried plotting datetime versus the along track distance, and nothing weird was happening there. Both tracks were recorded at the same time.
* Increased size of bounding polygon
* Extracted tracks next to 559 and 574
-----
* Made a notebook which extracts and saves selected tracks
  * Immensely improves processing speeds and prevents the kernel from regularly dying on me
* Successfully got cross correlation working in km.
* Visualized track proximity.
* Added a cropping mechanic which should mitigate the impact of track proximity?
* Added check on intersection distance accuracy!
* Noticed correlation between surface features and slope
* Identify points of high slope change
* Read a shitton
* Notice things correlate much better in areas with higher feature importances
* Created a reference track path, but im not sure how it would be usefull at all