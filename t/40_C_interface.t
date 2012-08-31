# tests for C interface
use strict;
use warnings;
use Test::More tests => 1;
use PDL;
#use Inline 'Noclean';
use Inline 'Pdlpp';
use PDL::Fit::ExpRate;

# Householder test
my $expected = zeroes(3)->grandom;
my $A = zeroes(3,3)->grandom;
my $y = ($A x $expected->transpose)->squeeze;
my $coefs1 = three_by_three_Householder($A, $y);
my $coefs2 = PDL::my_three_by_three_Householder($A, $y);
ok(all (approx($coefs1, $coefs2)), '3x3 Householder works')
	or diag("Expected $coefs1 but got $coefs2");




__DATA__

__Pdlpp__

use blib;
use PDL::PP::Include::Fit::ExpRate;

pp_def('my_three_by_three_Householder',
	Pars => 'A(n=3,m=3); y(n=3); [o] coef(n=3)',
	GenericTypes => ['D'],
	Code => pp_line_numbers(__LINE__, q{
		double tmp_A [3][3];
		double tmp_y [3];
		double tmp_c [3];
		threadloop %{
			/* Copy the values into the temporary matrices */
			loop (n) %{
				tmp_y[n] = $y();
				loop (m) %{
					tmp_A[n][m] = $A();
				%}
			%}
			
			/* Call the C-method */
			exprate_three_by_three_Householder(tmp_A, tmp_y, tmp_c);
			loop(n) %{
				$coef() = tmp_c[n];
			%}
		%}
	}),
);

