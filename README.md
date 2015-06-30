
# DNA_methylation
usage : "../sam2lollipop.pl --sam $sam --converted 0 or 1 --random 0 or 1";

The programs takes a sam file and search for the reads with at least one convertion (--converted 1) or at least one non-conversion (--converted 0) and plot the profiles of conversion for those reads (the first 100 of them). if --random 1 : a randomized profile based on the real data is ploted using the same criteria (--converted 1 or 0). 

