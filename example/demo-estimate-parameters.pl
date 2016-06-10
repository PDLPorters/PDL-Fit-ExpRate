# Illustrates a single Newton Method
use strict;
use warnings;

use PDL;
use PDL::Fit::ExpRate;

 # Make some noisy data
 my $xs = sequence(50);
 my $ys = 5 + 3 * exp($xs * 0.04) + grandom($xs) * 0.1;
 print "        A      B   lambda    err\n";
 
 # Print the original parameters
 printf "gave: %1.3f  %1.3f  %1.3f\n", 5, 3, 0.04;
 
 # Get estimates for parameters
 my ($A, $B, $lambda, $sum_sq_err) = exp_fit_estimate_parameters($xs, $ys);
 printf "est:  %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;
 
