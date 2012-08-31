# tests for quadratic fitting
use strict;
use warnings;
use Test::More tests => 11;
use PDL;
use PDL::Fit::ExpRate;

# Create a 3x3 diagonal matrix and make sure the we get the right answer
my $xs = sequence(1000);
my $ys = 2 + 4 * $xs + 6 * $xs**2;
my $coefs = fit_quadratic($xs, $ys);
my $expected = pdl(2, 4, 6);
ok(all (approx($coefs, $expected, 1e-2)), 'Super-simple quadratic fitting works')
	or diag("Expected $expected but got $coefs");

my $PGP_is_imported = 0;
# to-do: add a little noise and make sure it works
for (1..10) {
	use PDL::NiceSlice;
	$expected = zeroes(3)->grandom;
	$ys = $expected(0) + $expected(1) * $xs + $expected(2) * $xs**2
		+ $xs->grandom * $expected(1) / 10;
	$coefs = fit_quadratic($xs, $ys);
	my $const_diff = abs($coefs(0) - $expected(0));
	my $pct_diff = abs(($coefs(1:) - $expected(1:)) / $expected(1:));
	# Tolerance for constant is very loose because the noise makes it much less
	# robust. The others can be much tighter.
	ok( ( all($pct_diff->abs < pdl(0.001, 0.001)) and $const_diff < 0.1)
		, 'Quadratic fitting with noisy input works')
		or diag("Got $coefs but expected $expected\n"
			. "Got pct differences of $pct_diff and const diff of $const_diff");
}
