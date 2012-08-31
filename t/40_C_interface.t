# tests for C interface
use strict;
use warnings;
use Test::More tests => 10;
use PDL;
#use Inline 'Noclean';
use Inline 'Pdlpp';
use PDL::Fit::ExpRate;

# Householder (1)
my $expected = zeroes(3)->grandom;
my $A = zeroes(3,3)->grandom;
my $y = ($A x $expected->transpose)->squeeze;
my $coefs1 = three_by_three_Householder($A, $y);
my $coefs2 = PDL::my_three_by_three_Householder($A, $y);
ok(all (approx($coefs1, $coefs2)), '3x3 Householder works')
	or diag("Expected $coefs1 but got $coefs2");

# Quadratic fit (1)
my $xs = sequence(1000);
my $ys = 2 + 4 * $xs + 6 * $xs**2;
$coefs1 = fit_quadratic($xs, $ys);
$coefs2 = PDL::my_fit_quadratic($xs, $ys);
ok(all (approx($coefs1, $coefs2)), 'Quadratic fit works')
	or diag("Expected $coefs1 but got $coefs2");

# Exponential parameter estimate (4)
$xs = sequence(100);
$ys = 2 + 4 * exp($xs / -10);
my ($A1, $B1, $lambda1, $sum_sq_error1) = exp_fit_estimate_parameters($xs, $ys);
my ($A2, $B2, $lambda2, $sum_sq_error2)
	= PDL::my_exp_fit_estimate_parameters($xs, $ys);
is($A1, $A2, 'Estimates A the same');
is($B1, $B2, 'Estimates B the same');
is($lambda1, $lambda2, 'Estimates lambda the same');
is($sum_sq_error1, $sum_sq_error2, 'Calculates the same sum_sq_error');

# Newton step (4)
($A1, $B1, $lambda1, $sum_sq_error1) = exp_fit_newton_method_step(
	$xs, $ys, 0.1, $A1, $B1, $lambda1);
($sum_sq_error2) = PDL::my_exp_fit_newton_method_step(
	$xs, $ys, 0.1, $A2, $B2, $lambda2);
is($A1, $A2, 'Updates A the same');
is($B1, $B2, 'Updates B the same');
is($lambda1, $lambda2, 'Updates lambda the same');
is($sum_sq_error1, $sum_sq_error2, 'Calculates the same sum_sq_error');


__DATA__

__Pdlpp__

use blib;
use PDL::PP::Include::Fit::ExpRate 'PDL';

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

pp_def('my_fit_quadratic',
	Pars => 'xs(n); ys(n); [o] coefs(m=3)',
	GenericTypes => ['D'],
	Code => pp_line_numbers(__LINE__, q{
		_fit_quadratic($P(xs), $P(ys), $SIZE(n), $P(coefs));
	}),
);

pp_def('my_exp_fit_estimate_parameters' =>
	Pars => 'xs(n); ys(n); [o] As(); [o] Bs(); [o] lambdas(); [o] sum_sq_errors()',
	GenericTypes => ['D'],
	Code => q{
		/* Really simple: just call the C function */
		exprate_estimate_parameters($P(xs), $P(ys), $SIZE(n),
				$P(As), $P(Bs), $P(lambdas), $P(sum_sq_errors));
	},
);

pp_def('my_exp_fit_newton_method_step' =>
	Pars => 'xs(n); ys(n); trust_radii(); As(); Bs(); lambdas();
			[o] sum_sq_err()',
	GenericTypes => ['D'],
	Code => q{
		/* The A/B/lambda values are modified in-place */
		_exp_fit_newton_method_step($P(xs), $P(ys), $SIZE(n), $trust_radii(),
			&($As()), &($Bs()), &($lambdas()), &($sum_sq_err()));
	},
);
