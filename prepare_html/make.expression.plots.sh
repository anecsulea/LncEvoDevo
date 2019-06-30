#!/bin/bash

export sp=$1

############################################################################

export path=LncEvoDevo

export pathScripts=${path}/scripts/prepare_html

############################################################################
	
echo "#!/bin/bash" >  ${pathScripts}/bsub_script_exp_plot

echo "Rscript --vanilla ${pathScripts}/make.expression.plots.R ${sp}" >> ${pathScripts}/bsub_script_exp_plot

qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_exp_plots_${sp}.txt -e ${pathScripts}/std_error_exp_plots_${sp}.txt ${pathScripts}/bsub_script_exp_plot


############################################################################
