# Illustrates a single Newton Method
use strict;
use warnings;

use PDL;
use PDL::Fit::ExpRate;

 # Make some noisy data
 my $xs = sequence(500);
 my $ys = 5 + 3 * exp($xs * -0.05) + grandom($xs) * 0.5;
 print "   A      B      tau\n";
 
 # Use default parameters
 my ($A, $B, $tau) = fit_exp_rate($xs, $ys);
 printf " %1.3f  %1.3f  %1.3f  -> defaults\n", $A->sclr, $B->sclr, $tau->sclr;
 
 # See how things as we tighten the threshold a bit
 for my $threshold (0.1, 0.01, 0.001, 0.0001, 0.00001) {
     ($A, $B, $tau) = fit_exp_rate($xs, $ys,
       threshold => $threshold
     );
     printf " %1.3f  %1.3f  %1.3f  -> threshold = $threshold\n", $A->sclr, $B->sclr, $tau->sclr;
 }
 
 print '=' x 35, "\n";
 printf " %1.3f  %1.3f  %1.3f  -> values for generating data\n", 5, 3, -1/0.05;
 
 # Use callback to track the convergence of A
 my @As;
 ($A, $B, $tau) = fit_exp_rate($xs, $ys,
     run_each_iteration => sub {
         my %info = @_;
         push @As, $info{A};
     }
 );
 print "\nA:   ", join ("\n  -> ", @As), "\n";
 
 # Use callback to print diagnostics
 print '-' x 30, "\n";
 ($A, $B, $tau) = fit_exp_rate($xs, $ys,
     iterations => 100,
	 run_each_iteration => sub {
         my %info = @_;
         return if $info{round} % 5 != 0;
         print "Round $info{round}:\n";
         delete $info{round};
         for my $k (sort keys %info) {
             print "  $k: $info{$k}\n";
         }
     },
 );
 
 # Use callback to get more information at
 # the end of the fit
 print '-' x 30, "\n";
 ($A, $B, $tau) = fit_exp_rate($xs, $ys,
     iterations => 100,
     run_each_fit => sub {
         my %info = @_;
         print "Full details at end of fit:\n";
         for my $k (sort keys %info) {
             print " $k => $info{$k}\n";
         }
     }
 );
