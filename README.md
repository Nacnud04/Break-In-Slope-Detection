# Break-In-Slope-Detection

This repository adequately identifies the break in slope of an Ice Sheet-Shelf connection using ICESat-2 ATL06 data. <br>
A general outline of the repository and each notebook is contained here. More complex notebooks have detailed supporting documentation within the notebooks themselves.

## Simplified Process
* Compute flowslope
  * Performed by `S3-InitProcess.ipynb` using support files in the `ICESat2GroundingLineMigration` folder.
* Find break in slope & ice plains via flowslope
  * Computed by `FindSlopeBreak.ipynb`, `IbLine.ipynb` and `IbLine2.ipynb` via support files in the `src` folder.

