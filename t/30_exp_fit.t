# tests for quadratic fitting
use strict;
use warnings;
use Test::More tests => 2;
use PDL;
use PDL::Fit::ExpRate;

# Create some data and fit it
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
my $tau = -1e5;
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
