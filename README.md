# ToNER
ToNER is a tool for statistical modeling of enrichment from RNA-seq data comprising enriched and unenriched control libraries.

detail : http://www4a.biotec.or.th/GI/tools/toner (unavailable)

Requirements
•	python version 2.7 
    modules python
      -	matplotlib version 1.8.0 or higher
      -	numpy version 1.4.3 or higher
      -	scipy version 0.17.0 or higher
•	samtools version 1.2 or higher

Quick start 
>>> python ToNER.1.3.py –i folder_input [option]

Option program
-h, --help              	    Show help message and exit
-I, --input <string>	        Destination directory of input file ***requirement*** 	
                              The data must be already aligned to a reference genome sequence in BAM format. In 1 dataset must have 2 libraries, 
                              file name of enrich library must end with _1.bam and file name of control library must end with _2.bam such as ecoli_1.bam and 
                              ecoli_2.bam. If in input-folder have multiple datasets ToNER can analyze with meta-analysis for increase power of detection.
-o, --output	                Name directory of output file (default= ToNER_year_month_day_time)
-r, --reads <string>	        Position of reads for calculate [start , all , end] (default= all)
-p, --p_value <float>	        P-value cut off (default=0.05)
--none_pseudo <bool>	        None pseudo count 
-g, --gene <bool> 	          Add gene annotation information to analyzed results [start , end] requirement : gff annotation file in input folder
-t, –-total_read <int>	      filter reads count with total reads (sum of both libraries)
-e, --each_library <int>	    filter reads count with minimum read depth in either library
-q, --qqplot <float>	        cutoff r-squared to pass normal distribution after Box-cox transform (default=0.9)
-d, --distribution <string>	  type of distribution [normal, top_rank] (default= normal)
-c, --combine	                meta-analysis by Fisher’s combined probability test to calculate combined p-values
-s, --consensus	              meta-analysis by called significant among at least a minimum number of replicates specified
--less_memory <bool>	        In case of program has problem with memory, this option maybe can process but program must use long time (default=False)
--raw_data <bool>	            save depth and ratio to text file (default=True)
