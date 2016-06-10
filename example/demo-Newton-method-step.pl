# Illustrates a single Newton Method
use strict;
use warnings;

use PDL;
use PDL::Fit::ExpRate;

 # Make some noisy data
 my $xs = sequence(50);
 my $ys = 5 + 3 * exp($xs * 0.04) + grandom($xs) * 0.1;
 print "          A      B   lambda   err\n";
 
 # Get initial estimates for parameters
 my ($A, $B, $lambda, $sum_sq_err) = exp_fit_estimate_parameters($xs, $ys);
 printf "step 0: %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;
 
 # Perform many single steps to see how it goes
 my $trust_radius = 0.1;  # trust radius of 10%
 for (1 .. 9) {
   ($A, $B, $lambda, $sum_sq_err)
     = exp_fit_newton_method_step($xs, $ys, $trust_radius, $A, $B, $lambda);
   printf "step $_: %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;
 }
 
 # Print the original parameters
 print '=' x 35, "\n";
 printf "gave:   %1.3f  %1.3f  %1.3f\n", 5, 3, 0.04;
