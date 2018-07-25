# Priming of reversals and stability in perceptually coupled multistable displays

Data and analysis for the "Priming of reversals and stability in perceptually coupled multistable displays" manuscript and conference poster. Please refer to 
**Analysis.Rmd** for the analysis script and **Analysis.md** to view the results online.

# Data format
* `Observer` A unique participant ID 
* `SessionID` Session timestampt
* `Block` Overall block index
* `Trial` Trial index within the block
* `OnsetDelay` Randomized display onset delay in seconds. 
* `Yellow` Logical, `True` means that **left** object was colored yellow and the right object was white, `False` means that the left one was white and the right one was yellow.
* `Switch` Type of the exogenous trigger. Can be either `neither` (objects' motion remained unperturbed), `left`, `right` (trigger was applied to a single object), or `both` (trigger was applied to both objects).
* `SwitchTime` Time in seconds relative to the display onset when the trigger was applied. Irrelevant for `neither` trigger type.
* `Gap` Half-distance between objects in the units of their width. Negative gap means _overlap_, zero gap - _touching_, positive value - _gap_ layout.
* `Direction` Formal direction of rotation of the ideal sphere (irrelevant, as both displays were fully ambiguous).
* `Response` Participant's response on perceived stability of both objects. Can be `neither` (both spheres remained stable), `left`, `right` (a single sphere reversed its direction of rotation), or `both` (both spheres reversed the rotation).
* `RT` Response time in seconds.

## License
All data (and associated content) is licensed under the [CC-By Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/). All code is licensed
under the [MIT License](http://www.opensource.org/licenses/mit-license.php).
