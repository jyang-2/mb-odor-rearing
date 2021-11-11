
This folder (step00_process_suite2p) contains the outputs of running step00_process_suite2p.py on the initial suite2p outputs.

It should contain the following files:
	- params.json: contains the parameters used for extracting C_dec, the cells used, and the cost
			- dcnv_params
			- oasis_params
			- respvec_ops 
			- good_cells: cells used to compute cost function, specified by user
	
	- model_selection.json: contains a list of "models" (or parameter combos) to try, and their cost. 
							The selected parameters in params.json come from the "model" with the lowest cost.
	
	- Fcc: outputs from baseline correction w/ suite2p.dcnv, run block by block (using df_frametimeinfo) before np.hstack
 - Fq: quantile-normed Fcc (preprocess.RobustScaler)
 - C_dec: deconvolved traces
 - Sp : Spike array	- oasis_results: final output of deconvolution
	-    
