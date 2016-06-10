# Illustrates a single Newton Method
use strict;
use warnings;

use PDL;
use PDL::Fit::ExpRate;

 # Make some noisy data
 my $xs = sequence(500);
 my $ys = 5 + 3 * exp($xs * -0.05) + grandom($xs) * 0.5;
 print "          A      B   lambda   err\n";
 
 # Get initial estimates for parameters
 my ($A, $B, $lambda, $sum_sq_err) = exp_fit_estimate_parameters($xs, $ys);
 printf "step 0: %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;
 
 # Perform many single steps to see how it goes
 my $trust_radii = pdl(0.01, 0.03, 0.1, 0.3, 1, 3);
 
 ($A, $B, $lambda, $sum_sq_err)
   = exp_fit_newton_method_step($xs, $ys, $trust_radii, $A, $B, $lambda);
 print cat($trust_radii, $A, $B, $lambda, $sum_sq_err)->transpose;
 
 # Print the original parameters
 print '=' x 35, "\n";
 printf "gave:   %1.3f  %1.3f  %1.3f\n", 5, 3, 0.04;
