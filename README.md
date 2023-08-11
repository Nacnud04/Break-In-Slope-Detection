# Break-In-Slope-Detection

This repository adequately identifies the break in slope of an Ice Sheet-Shelf connection using ICESat-2 ATL06 data. <br>
A general outline of the repository and each notebook is contained here. More complex notebooks have detailed supporting documentation within the notebooks themselves.

## Simplified Process
* Compute flowslope
  * Performed by `S3-InitProcess.ipynb` using source files in the `ICESat2GroundingLineMigration` folder.
* Find break in slope & ice plains via flowslope
  * Computed by `FindSlopeBreak.ipynb`, `IbLine.ipynb` and `IbLine2.ipynb` via source files in the `src` folder.

## Main Directory
### Folders
Each folders purpose is decribed by `FILESTRUCT.md` in the `Documentation` folder.

### Notebooks
Each notebook's function is described below:

**FindSlopeBreak**: Imports a single beam of a cycle granule & cycle from the `Saves` folder and computes the predicted point of slope break. This is highly useful for debugging or understanding why the algorithm operates as it does, as all possible picks are reported, as well as the parameters which ended up resulting in the final decision. If the algorithm is not confident that there is a break in slope it will error.

**GenDirs**: When data is produced by `S3-InitProcess` it needs to be exported into pre-existing directories. This notebook produces those directories for a given region across a predefined range of cycles.

**IbCombined**: When break in slope predictions are exported by `IbLine2` they are output into the `Line` directory. This notebook imports all exported data and displays it all in a variety of ways. This includes: Comparisons with other grounding line models, visualizing the quality score parameter, median location of a single break in slope across cycles, as well as displaying the raw data.

**IbLine2**: Relies on `src/SlopeBreak2.py` (similarly `IbLine` relies on `src/SlopeBreak.py`). This notebook imports an entire cycle of data exported by `S3-InitProcess.py` and finds the break in slope for each beam and ground track. Then it exports the data to `Line`

