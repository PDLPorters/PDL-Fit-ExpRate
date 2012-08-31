# tests for quadratic fitting
use strict;
use warnings;
use Test::More tests => 2;
use PDL;
use PDL::Fit::ExpRate;

# Create a 3x3 diagonal matrix and make sure the we get the right answer
my $xs = sequence(30);
my $ys = 150 + 10 * exp($xs / -10);
my ($As, $Bs, $taus) = fit_exp_rate($xs, $ys
	, threshold => 0.00001
	, iterations => 200
	, trust_radius => 0.1,
);
my $expected = pdl(150, 10, -10);
my $coefs = pdl($As, $Bs, $taus)->flat;
ok(all (approx($coefs, $expected, 1e-2)), 'Super-simple exponential fitting works')
	or diag("Got $coefs but expected $expected");

$xs = sequence(30000);
my $tau = -1e7;
$ys = 150 + 10 * exp($xs / $tau);
($As, $Bs, $taus) = fit_exp_rate($xs, $ys
	, threshold => 0.0001
	, iterations => 400
	, trust_radius => 0.1,
);
$coefs = pdl($As, $Bs, $taus)->flat;
$expected = pdl(150, 10, $tau);
ok(all (abs(($coefs - $expected) / $expected) < 1e-7)
	, 'Super-simple exponential fitting with lots of data works')
	or diag("Got $coefs but expected $expected");


__END__
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
