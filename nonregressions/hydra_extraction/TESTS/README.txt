
This folder contains a series of non-regression test for XRS data extraction from ID20.

All tests are based on the same data set taken during "run3_16/run4_es431" of a small stishovite SiO2 (reference sample for an official experiment).

test1: auto_ROIs_pixel.yaml
	Automatic selection of ROIs, pixel-by-pixel compensation.

test2: linear_ROIs_pixel.yaml
	Line ROIs selection, pixel-by-pixel compensation.

test3: polygon_ROIs_pixel.yaml
	Chosing ROIs using polygons, pixel-by-pixel compensation.

test4: zoom_ROIs_NNMF_pixel.yaml
	Zoom ROIs selection, refinement of the ROIs using non-negative matrix factorization based on a single oxygen K-edge scan, pixel-by-pixel 	compensation.

test5: zoom_ROIs_pixel.yaml
	Zoom ROI selection, pixel-by-pixel compensation.

test6: zoom_ROIs_sum.yaml
	Zoom ROI selection, simple sum of signals within each ROI, no compensation.

test7: hydra_extraction.yaml
	Generic example.




